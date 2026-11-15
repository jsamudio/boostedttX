#!/usr/bin/env python3
"""
PBS Pro submission tools for the pocket-coffea manual-jobs workflow (see executors_pbs.py).

Pure Python (no pocket-coffea import): runs on the host, where qsub is available, as well
as inside the apptainer image.

A jobs directory <outputdir>/job is prepared by
    pocket-coffea run --cfg <config> -o <outputdir> --executor pbs \
        --executor-custom-setup ../common/executors_pbs.py -ro run_options_pbs.yaml
and contains one pickled configurator per job (config_job_N.pkl) and its fileset (fileset_job_N.json),
job.sh, jobs_config.yaml and pbs_settings.json.
This module turns it into a PBS job array (one subjob per pocket-coffea job, each running
the futures executor on all cores of its node inside apptainer) plus a merge that runs after
the array.

Merge modes ('merge-mode' in run_options_pbs.yaml / pbs_settings.json):
  dataset (default)  one merge per dataset, as a PBS array (merge_group.py), writing
                     <outputdir>/merged/<dataset>.coffea: each a complete pocket-coffea output of its
                     dataset (columns, genweight and theory sums, cutflow), so a merge needs the memory
                     of one dataset, not of the whole group. Jobs sharing a dataset are merged together
                     (merge_groups). merged_outputs(<outputdir>) returns the files to read.
  all                one merge of every job into <outputdir>/output_all.coffea (the old behaviour).

    python pbs_tools.py submit   <outputdir>    # array + merge
    python pbs_tools.py pending  <outputdir>... # submit the prepared, not yet submitted dirs, until
                                                # the per-user queue limit is reached (rerun later)
    python pbs_tools.py status   <outputdir>    # done / failed / running / idle per job, merged datasets
    python pbs_tools.py resubmit <outputdir>    # failed jobs (+ --include-unfinished), then the merges
                                                # of their datasets
    python pbs_tools.py merge    <outputdir>    # (re)submit the merges whose output is missing
                                                # (--all: every dataset again)
    python pbs_tools.py groups   <outputdir>    # store the dataset groups in pbs_settings.json; needs
                                                # PyYAML: for dirs prepared before the per-dataset merge,
                                                # run it once inside the image if the host has no PyYAML
"""
import argparse
import os
import shutil
import subprocess
import json
import sys

SETTINGS_FILE = "pbs_settings.json"
SUBMISSIONS_FILE = "pbs_submissions.json"
# Set by run_pipeline.py: the executor only prepares the jobs, submission happens on the host
NO_SUBMIT_ENV = "PC_PBS_NO_SUBMIT"

# Defaults, overridable in run_options_pbs.yaml
DEFAULTS = {
    "queue": "hep",
    "cores-per-worker": 8,
    "mem-per-worker": "32gb",          # whole job (all cores)
    "walltime": None,                  # e.g. "12:00:00"; None = no request
    "merge-mem": "32gb",               # per merge job (per dataset in merge-mode 'dataset')
    "merge-mode": "dataset",           # 'dataset': merged/<dataset>.coffea; 'all': output_all.coffea
    "scratch-dir": "/tmp",             # job-local working directory (must be bound in the image)
    "apptainer-image": "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-latest",
    "apptainer-binds": ["/cvmfs", "/cms/data", "/home/$USER", "/tmp"],
}


# ---------------------------------------------------------------- helpers
def find_jobs_dir(path):
    """The prepared jobs directory of an output dir, or None if it was not prepared."""
    for d in (path, os.path.join(path, "job")):
        if os.path.isfile(os.path.join(d, SETTINGS_FILE)):
            return d
    return None


def jobs_dir_of(path):
    """Accept either the pocket-coffea output dir or its job/ subdirectory."""
    for d in (path, os.path.join(path, "job")):
        if os.path.isfile(os.path.join(d, SETTINGS_FILE)):
            return os.path.abspath(d)
    sys.exit(f"No {SETTINGS_FILE} found in {path} or {path}/job")


def load_settings(jobs_dir):
    with open(os.path.join(jobs_dir, SETTINGS_FILE)) as f:
        return json.load(f)


def apptainer_prefix(settings, env=""):
    """`apptainer exec` with the binds and KRB5CCNAME; `env` = extra `--env` options (before the image)."""
    binds = " ".join(f"-B {b}" for b in settings["apptainer-binds"])
    return f'apptainer exec {binds} --env KRB5CCNAME="${{KRB5CCNAME:-}}" {env + " " if env else ""}{settings["apptainer-image"]}'


def _resources(settings, ncpus, mem):
    lines = [f"#PBS -q {settings['queue']}", f"#PBS -l select=1:ncpus={ncpus}:mem={mem}"]
    if settings.get("walltime"):
        lines.append(f"#PBS -l walltime={settings['walltime']}")
    return "\n".join(lines)


def _mem_gb(mem):
    """'32gb' / '32GB' / '32000mb' -> GB (float)."""
    m = str(mem).strip().lower()
    for unit, scale in (("tb", 1024.), ("gb", 1.), ("mb", 1 / 1024.), ("kb", 1 / 1024. ** 2)):
        if m.endswith(unit):
            return float(m[:-len(unit)]) * scale
    return float(m) / 1024. ** 3        # plain bytes


class QueueFull(Exception):
    """qsub refused the job: the per-user limit of queued jobs is reached."""


def qsub(script, depend=None):
    if shutil.which("qsub") is None:
        sys.exit("qsub not found. Inside the apptainer image PBS is not available: run\n"
                 f"    python {os.path.abspath(__file__)} submit <outputdir>\n"
                 "on the host (outside the container).")
    cmd = ["qsub"] + (["-W", f"depend={depend}"] if depend else []) + [script]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        if "would exceed" in res.stderr:
            raise QueueFull(res.stderr.strip())
        sys.exit(f"qsub failed ({' '.join(cmd)}):\n{res.stderr}")
    return res.stdout.strip()


def record_submission(jobs_dir, entry):
    path = os.path.join(jobs_dir, SUBMISSIONS_FILE)
    subs = []
    if os.path.exists(path):
        with open(path) as f:
            subs = json.load(f)
    subs.append(entry)
    with open(path, "w") as f:
        json.dump(subs, f, indent=1)


def job_ids(jobs_dir):
    """pocket-coffea job ids job_0 ... job_<N-1> (one chunksize per job in the settings)."""
    return list(range(len(load_settings(jobs_dir)["chunksizes"])))


# ---------------------------------------------------------------- per-dataset merge
MERGED_DIR = "merged"            # <outputdir>/merged/<dataset>.coffea (also in merge_group.py)


def merge_mode(settings):
    return settings.get("merge-mode") or DEFAULTS["merge-mode"]


def group_jobs(job_datasets):
    """Groups of jobs to merge together: {name: [job ids]}, one per dataset. job_datasets is
    {job id: [datasets of the job]}; jobs that share a dataset end up in the same group (connected
    components), so every merged file holds complete datasets. Name: the dataset, or the first
    dataset + '_plus<N>' if a group holds several."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, datasets in job_datasets.items():
        if not datasets:
            raise ValueError(f"job {i} has no datasets")
        for d in datasets[1:]:
            parent[find(d)] = find(datasets[0])
        find(datasets[0])
    comps = {}
    for i, datasets in job_datasets.items():
        comps.setdefault(find(datasets[0]), set()).add(i)
    groups = {}
    for jobs in comps.values():
        datasets = sorted({d for i in jobs for d in job_datasets[i]})
        name = datasets[0] + (f"_plus{len(datasets) - 1}" if len(datasets) > 1 else "")
        name = "".join(c if c.isalnum() or c in "_-." else "_" for c in name)
        groups[name] = sorted(int(i) for i in jobs)
    return dict(sorted(groups.items()))


def merge_groups(jobs_dir):
    """{group name: [job ids]}: from pbs_settings.json ('merge-groups', written at preparation by
    executors_pbs.py), else derived from jobs_config.yaml (needs PyYAML) and stored there."""
    settings = load_settings(jobs_dir)
    if settings.get("merge-groups"):
        return settings["merge-groups"]
    try:
        import yaml
    except ImportError:
        sys.exit(f"{jobs_dir}: no dataset groups in {SETTINGS_FILE} and no PyYAML here to read "
                 f"jobs_config.yaml. Run once inside the apptainer image:\n"
                 f"    python {os.path.abspath(__file__)} groups {os.path.dirname(jobs_dir)}")
    with open(os.path.join(jobs_dir, "jobs_config.yaml")) as f:
        jobs_list = yaml.safe_load(f)["jobs_list"]
    groups = group_jobs({int(name.split("_")[1]): sorted((job.get("filesets") or {}).keys())
                         for name, job in jobs_list.items()})
    settings["merge-groups"] = groups
    with open(os.path.join(jobs_dir, SETTINGS_FILE), "w") as f:
        json.dump(settings, f, indent=1)
    return groups


def merged_file(settings, group):
    return os.path.join(settings["output-dir"], MERGED_DIR, f"{group}.coffea")


def merged_outputs(path, allow_incomplete=False):
    """The merged .coffea files of an output dir, to pass to the analysis scripts:
    merged/<dataset>.coffea (per-dataset merge) or output_all.coffea (single merge, or a local run).
    Raises if both kinds are present (the same events twice) or, for a PBS dir, if datasets are
    still unmerged (unless allow_incomplete). Returns [] if nothing is merged yet."""
    path = os.path.abspath(path)
    merged_dir = os.path.join(path, MERGED_DIR)
    per_ds = (sorted(os.path.join(merged_dir, f) for f in os.listdir(merged_dir) if f.endswith(".coffea"))
              if os.path.isdir(merged_dir) else [])
    single = os.path.join(path, "output_all.coffea")
    if per_ds and os.path.exists(single):
        raise RuntimeError(f"{path} has both output_all.coffea and {MERGED_DIR}/*.coffea (the same events "
                           f"twice): remove one of them")
    if not per_ds:
        return [single] if os.path.exists(single) else []
    jobs_dir = find_jobs_dir(path)
    if jobs_dir is not None and not allow_incomplete:
        settings = load_settings(jobs_dir)
        groups = settings.get("merge-groups")
        if groups:
            # by file name inside merged/ (the stored output-dir may differ from `path` by symlinks)
            names = {os.path.basename(f) for f in per_ds}
            missing = [g for g in groups if f"{g}.coffea" not in names]
            if missing:
                raise RuntimeError(f"{path}: {len(missing)} of {len(groups)} datasets not merged yet "
                                   f"(e.g. {missing[:3]}); see pbs_tools.py status / merge")
            extra = sorted(names - {f"{g}.coffea" for g in groups})
            if extra:
                raise RuntimeError(f"{path}: files in {MERGED_DIR}/ that are not a dataset group of this "
                                   f"dir: {extra}")
    return per_ds


# ---------------------------------------------------------------- scripts
def write_joblist(jobs_dir, name, ids, settings):
    """One line per array index: '<pocket-coffea job id> <chunksize>'."""
    chunks = settings["chunksizes"]
    path = os.path.join(jobs_dir, name)
    with open(path, "w") as f:
        for i in ids:
            f.write(f"{i} {chunks[i]}\n")
    return path


def write_array_script(jobs_dir, settings, joblist, tag):
    n = sum(1 for _ in open(joblist))
    array_line = f"#PBS -J 0-{n - 1}" if n > 1 else ""
    script = f"""#!/bin/bash
#PBS -N {settings['name']}
{_resources(settings, settings['cores-per-worker'], settings['mem-per-worker'])}
{array_line}
#PBS -j oe
#PBS -o {jobs_dir}/logs/

# One pocket-coffea job per array index: futures executor on all cores of the node.
JOBS_DIR={jobs_dir}
IDX=${{PBS_ARRAY_INDEX:-0}}
read JOBID CHUNKSIZE <<< "$(sed -n "$((IDX + 1))p" {joblist})"

WORK={settings['scratch-dir']}/$USER/pc_{settings['name']}_${{JOBID}}_$$
mkdir -p "$WORK" && cd "$WORK" || exit 1
# pocket-coffea >= 0.9.13 job.sh: runs `--cfg {SHARED_CONFIG_NAME} --fileset "$2"`; our per-job
# configurator (restricted to the job's samples) is copied under that name
cp "$JOBS_DIR/config_job_$JOBID.pkl" {SHARED_CONFIG_NAME} || exit 1
cp "$JOBS_DIR/fileset_job_$JOBID.json" "$JOBS_DIR/job.sh" "$JOBS_DIR/{settings['inner-yaml']}" . || exit 1

{apptainer_prefix(settings, '--env _CONDOR_SCRATCH_DIR="$WORK"')} \\
    bash job.sh "$JOBID" "fileset_job_$JOBID.json" "$CHUNKSIZE" {settings['cores-per-worker']} > "$JOBS_DIR/logs/job_$JOBID.log" 2>&1
STATUS=$?

cd / && rm -rf "$WORK"
exit $STATUS
"""
    path = os.path.join(jobs_dir, f"pbs_array_{tag}.sh")
    with open(path, "w") as f:
        f.write(script)
    return path


def write_merge_script(jobs_dir, settings, tag, groups=None):
    """Merge job script: per-dataset array over `groups` ({name: [job ids]}) in merge-mode 'dataset',
    else the single merge into output_all.coffea."""
    mem_gb = _mem_gb(settings["merge-mem"])
    if merge_mode(settings) == "dataset":
        if not groups:
            raise ValueError("no dataset groups to merge")
        groups_file = os.path.join(jobs_dir, f"merge_groups_{tag}.txt")
        with open(groups_file, "w") as f:
            for name, ids in groups.items():
                f.write(f"{name} {','.join(str(i) for i in ids)}\n")
        array_line = f"#PBS -J 0-{len(groups) - 1}" if len(groups) > 1 else ""
        common = settings.get("common-dir", os.path.dirname(os.path.abspath(__file__)))
        script = f"""#!/bin/bash
#PBS -N {settings['name'][:10]}_merge
{_resources(settings, 1, settings['merge-mem'])}
{array_line}
#PBS -j oe
#PBS -o {jobs_dir}/logs/

# One merge per dataset group ({len(groups)}): merged/<group>.coffea; jobs with missing outputs are marked .failed
IDX=${{PBS_ARRAY_INDEX:-0}}
NAME=$(sed -n "$((IDX + 1))p" {groups_file} | cut -d' ' -f1)
cd {settings['analysis-dir']}
{apptainer_prefix(settings)} bash -c '
    export PYTHONPATH={settings['analysis-dir']}:{common}:$PYTHONPATH
    python {common}/merge_group.py {jobs_dir} {groups_file} '"$IDX"' -m {max(1.0, mem_gb / 2):.0f}
' > {jobs_dir}/logs/merge_{tag}_"$NAME".log 2>&1
"""
        path = os.path.join(jobs_dir, f"pbs_merge_{tag}.sh")
        with open(path, "w") as f:
            f.write(script)
        return path

    script = f"""#!/bin/bash
#PBS -N {settings['name'][:10]}_merge
{_resources(settings, 1, settings['merge-mem'])}
#PBS -j oe
#PBS -o {jobs_dir}/logs/

# Merge (and postprocess) all job outputs; jobs with missing outputs are marked .failed
cd {settings['analysis-dir']}
{apptainer_prefix(settings)} bash -c '
    export PYTHONPATH={settings['analysis-dir']}:{settings.get('common-dir', settings['analysis-dir'])}:$PYTHONPATH
    merge-outputs -jc {jobs_dir}/jobs_config.yaml -o {settings['output-dir']}/output_all.coffea \\
        --mark-failed -f -m {max(1.0, mem_gb / 2):.0f}
' > {jobs_dir}/logs/merge_{tag}.log 2>&1
"""
    path = os.path.join(jobs_dir, f"pbs_merge_{tag}.sh")
    with open(path, "w") as f:
        f.write(script)
    return path


# ---------------------------------------------------------------- grid proxy
PROXY_MIN_HOURS = 1        # refuse to submit below this
# Names used by the pocket-coffea (>= 0.9.13) job.sh: shared configurator in the work dir, and the
# per-job fileset passed as $2 (executors_manual_jobs.SHARED_CONFIG_NAME / fileset_job_N.json)
SHARED_CONFIG_NAME = "configurator_jobs.pkl"
FILESET_JOB = "fileset_job_{i}.json"
PROXY_WARN_HOURS = 24      # warn below this (jobs may queue / run long)


def refresh_proxy(settings):
    """Copy the current VOMS proxy to the path the jobs use (X509_USER_PROXY in job.sh) and check
    its lifetime, so that (re)submissions always ship a valid proxy. No-op without a proxy."""
    dest = settings.get("x509-path")
    if not dest:
        return
    src = os.environ.get("X509_USER_PROXY") or f"/tmp/x509up_u{os.getuid()}"
    if os.path.isfile(src) and os.path.abspath(src) != os.path.abspath(dest):
        shutil.copyfile(src, dest)
        os.chmod(dest, 0o600)
    if not os.path.isfile(dest):
        sys.exit(f"No VOMS proxy found ({src}). Run `voms-proxy-init -voms cms -rfc --valid 168:0` "
                 f"(e.g. inside the apptainer image) and submit again.")
    if shutil.which("voms-proxy-info") is None:
        print(f"  proxy: {dest} (voms-proxy-info not available here, lifetime not checked)")
        return
    res = subprocess.run(["voms-proxy-info", "-file", dest, "-timeleft"], capture_output=True, text=True)
    hours = int(res.stdout.strip() or 0) / 3600. if res.returncode == 0 else 0.
    if hours < PROXY_MIN_HOURS:
        sys.exit(f"VOMS proxy {dest} expired or invalid ({hours:.1f} h left). "
                 f"Run `voms-proxy-init -voms cms -rfc --valid 168:0` and submit again.")
    if hours < PROXY_WARN_HOURS:
        print(f"  [WARNING] VOMS proxy has only {hours:.1f} h left; jobs starting later will fail.")
    else:
        print(f"  proxy: {dest} ({hours:.0f} h left)")


# ---------------------------------------------------------------- commands
def submit(path, ids=None, tag="all"):
    """Submit the job array (all jobs, or `ids`) and a merge job that runs after it."""
    jobs_dir = jobs_dir_of(path)
    settings = load_settings(jobs_dir)
    refresh_proxy(settings)
    os.makedirs(os.path.join(jobs_dir, "logs"), exist_ok=True)
    ids = job_ids(jobs_dir) if ids is None else ids
    groups = groups_of(jobs_dir, settings, ids)       # resolve first: fail before touching job states or qsub
    for i in ids:
        for state in ("failed", "running", "done"):
            p = os.path.join(jobs_dir, f"job_{i}.{state}")
            if os.path.exists(p):
                os.remove(p)
        open(os.path.join(jobs_dir, f"job_{i}.idle"), "w").close()

    joblist = write_joblist(jobs_dir, f"joblist_{tag}.txt", ids, settings)
    array_id = qsub(write_array_script(jobs_dir, settings, joblist, tag))
    try:
        merge_id = qsub(write_merge_script(jobs_dir, settings, tag, groups), depend=f"afterany:{array_id}")
    except QueueFull:
        # the array is in; the merge job is submitted later by submit_pending
        record_submission(jobs_dir, {"tag": tag, "array": array_id, "merge": None, "jobs": list(ids)})
        print(f"{settings['name']}: submitted {len(ids)} job(s) as {array_id}; merge job not yet (queue full)")
        raise
    record_submission(jobs_dir, {"tag": tag, "array": array_id, "merge": merge_id, "jobs": list(ids)})
    target = (f"{len(groups)} dataset(s) -> {settings['output-dir']}/{MERGED_DIR}/" if groups
              else f"{settings['output-dir']}/output_all.coffea")
    print(f"{settings['name']}: submitted {len(ids)} job(s) as {array_id}; merge job {merge_id}, {target}")
    return array_id, merge_id


def groups_of(jobs_dir, settings, ids=None):
    """Dataset groups to merge after the jobs `ids` (all if None): those containing one of them, plus
    every group without a merged file yet (e.g. a dir merged before into output_all.coffea, or a merge
    that failed). None in merge-mode 'all'."""
    if merge_mode(settings) != "dataset":
        return None
    groups = merge_groups(jobs_dir)
    if ids is None:
        return groups
    ids = set(int(i) for i in ids)
    return {g: j for g, j in groups.items() if ids & set(j) or not os.path.exists(merged_file(settings, g))}


def submissions(jobs_dir):
    path = os.path.join(jobs_dir, SUBMISSIONS_FILE)
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def submit_pending(paths):
    """Submit, in order, every prepared output dir that has no job array submitted yet, and the
    missing merge jobs, until PBS refuses (per-user queue limit). Returns the number still pending."""
    pending = []
    for path in paths:
        jobs_dir = find_jobs_dir(path)
        if jobs_dir is None:
            continue
        arrays = [s for s in submissions(jobs_dir) if s.get("array")]
        if not arrays:
            pending.append((jobs_dir, "array"))
        elif arrays[-1].get("merge") is None:
            pending.append((jobs_dir, "merge"))
    for n, (jobs_dir, what) in enumerate(pending):
        try:
            if what == "array":
                submit(jobs_dir)
            else:
                submit_missing_merge(jobs_dir)
        except QueueFull as err:
            print(f"Queue limit reached ({err}).\n{len(pending) - n} submission(s) still pending: "
                  f"run the same command again later.")
            return len(pending) - n
    print("Nothing pending." if not pending else "All pending submissions done.")
    return 0


def submit_missing_merge(jobs_dir):
    """Merge job of the last array, which could not be submitted when the queue was full."""
    settings = load_settings(jobs_dir)
    subs = submissions(jobs_dir)
    last = [s for s in subs if s.get("array")][-1]
    active = subprocess.run(["qstat", last["array"]], capture_output=True).returncode == 0
    merge_id = qsub(write_merge_script(jobs_dir, settings, last["tag"], groups_of(jobs_dir, settings, last.get("jobs"))),
                    depend=f"afterany:{last['array']}" if active else None)
    last["merge"] = merge_id
    with open(os.path.join(jobs_dir, SUBMISSIONS_FILE), "w") as f:
        json.dump(subs, f, indent=1)
    print(f"{settings['name']}: merge job {merge_id}")


def states(jobs_dir):
    out = {}
    for i in job_ids(jobs_dir):
        state = "missing"
        for s in ("done", "failed", "running", "idle"):
            if os.path.exists(os.path.join(jobs_dir, f"job_{i}.{s}")):
                state = s
                break
        out[i] = state
    return out


def status(path):
    jobs_dir = jobs_dir_of(path)
    settings = load_settings(jobs_dir)
    st = states(jobs_dir)
    counts = {s: sum(1 for v in st.values() if v == s) for s in ("done", "failed", "running", "idle", "missing")}
    if merge_mode(settings) == "dataset" and settings.get("merge-groups"):
        groups = settings["merge-groups"]
        n = sum(os.path.exists(merged_file(settings, g)) for g in groups)
        merged = f"{n}/{len(groups)} datasets"
    elif merge_mode(settings) == "dataset":
        n = len([f for f in os.listdir(os.path.join(settings["output-dir"], MERGED_DIR)) if f.endswith(".coffea")]) \
            if os.path.isdir(os.path.join(settings["output-dir"], MERGED_DIR)) else 0
        merged = f"{n} datasets (groups not stored yet)"
    else:
        merged = "yes" if os.path.exists(os.path.join(settings["output-dir"], "output_all.coffea")) else "no"
    if os.path.exists(os.path.join(settings["output-dir"], "output_all.coffea")) and merge_mode(settings) == "dataset":
        merged += "; output_all.coffea from an older single merge is also present (remove one)"
    print(f"{settings['name']:24} " + "  ".join(f"{k}={v}" for k, v in counts.items()) + f"  merged={merged}")
    failed = [i for i, s in st.items() if s == "failed"]
    if failed:
        print(f"    failed jobs: {failed}  (logs: {jobs_dir}/logs/job_<id>.log)")
    return st


def _active(jobs_dir):
    """Array jobs of this directory still known to PBS as queued/running."""
    path = os.path.join(jobs_dir, SUBMISSIONS_FILE)
    if not os.path.exists(path) or shutil.which("qstat") is None:
        return []
    with open(path) as f:
        subs = json.load(f)
    return [s["array"] for s in subs if "array" in s
            if subprocess.run(["qstat", s["array"]], capture_output=True).returncode == 0]


def resubmit(path, include_unfinished=False):
    jobs_dir = jobs_dir_of(path)
    st = states(jobs_dir)
    redo = [i for i, s in st.items() if s == "failed"]
    if include_unfinished:
        active = _active(jobs_dir)
        if active:
            sys.exit(f"PBS still has active arrays {active} for this directory; wait for them "
                     f"(or qdel them) before resubmitting unfinished jobs.")
        redo += [i for i, s in st.items() if s in ("running", "idle", "missing")]
    if not redo:
        print("Nothing to resubmit.")
        return
    n = 1
    while os.path.exists(os.path.join(jobs_dir, f"joblist_r{n}.txt")):
        n += 1
    print(f"Resubmitting jobs {sorted(redo)}")
    submit(jobs_dir, ids=sorted(redo), tag=f"r{n}")


def merge(path, all_groups=False):
    """Submit the merge now (no dependency): in merge-mode 'dataset' the datasets whose merged file is
    missing (all of them with all_groups), else the single merge."""
    jobs_dir = jobs_dir_of(path)
    settings = load_settings(jobs_dir)
    groups = groups_of(jobs_dir, settings)
    if groups is not None and not all_groups:
        groups = {g: j for g, j in groups.items() if not os.path.exists(merged_file(settings, g))}
        if not groups:
            print(f"{settings['name']}: every dataset is merged (--all to merge them again)")
            return
    n = 1
    while os.path.exists(os.path.join(jobs_dir, f"pbs_merge_m{n}.sh")):
        n += 1
    merge_id = qsub(write_merge_script(jobs_dir, settings, f"m{n}", groups))
    record_submission(jobs_dir, {"tag": f"m{n}", "merge": merge_id})
    print(f"{settings['name']}: merge job {merge_id}"
          + (f" ({len(groups)} dataset(s): {', '.join(list(groups)[:5])}{' ...' if len(groups) > 5 else ''})"
             if groups else ""))


def _run_command(args):
    for d in args.outputdirs:
        if args.command == "submit":
            submit(d)
        elif args.command == "status":
            status(d)
        elif args.command == "resubmit":
            resubmit(d, args.include_unfinished)
        elif args.command == "groups":
            groups = merge_groups(jobs_dir_of(d))
            print(f"{d}: {len(groups)} dataset groups stored ({sum(len(j) for j in groups.values())} jobs)")
        else:
            merge(d, args.all)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["submit", "pending", "status", "resubmit", "merge", "groups"])
    ap.add_argument("outputdirs", nargs="+", help="pocket-coffea output dir(s) (or their job/ subdir)")
    ap.add_argument("--include-unfinished", action="store_true",
                    help="resubmit: also jobs left running/idle/missing (e.g. killed by the batch system)")
    ap.add_argument("--all", action="store_true",
                    help="merge: merge every dataset again, not only those without a merged file")
    args = ap.parse_args()
    if args.command == "pending":
        sys.exit(1 if submit_pending(args.outputdirs) else 0)
    try:
        _run_command(args)
    except QueueFull as err:
        sys.exit(f"Queue limit reached ({err}); try again later.")
