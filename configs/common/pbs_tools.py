#!/usr/bin/env python3
"""
PBS Pro submission tools for the pocket-coffea manual-jobs workflow (see executors_pbs.py).

Pure Python (no pocket-coffea import): runs on the host, where qsub is available, as well
as inside the apptainer image.

A jobs directory <outputdir>/job is prepared by
    pocket-coffea run --cfg <config> -o <outputdir> --executor pbs \
        --executor-custom-setup ../common/executors_pbs.py -ro run_options_pbs.yaml
and contains one pickled configurator per job, job.sh, jobs_config.yaml and pbs_settings.json.
This module turns it into a PBS job array (one subjob per pocket-coffea job, each running
the futures executor on all cores of its node inside apptainer) plus a merge job that runs
after the array and writes <outputdir>/output_all.coffea.

    python pbs_tools.py submit   <outputdir>    # array + merge job
    python pbs_tools.py status   <outputdir>    # done / failed / running / idle per job
    python pbs_tools.py resubmit <outputdir>    # failed jobs (+ --include-unfinished) and a new merge
    python pbs_tools.py merge    <outputdir>    # (re)submit only the merge job
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
    "queue": "batch",
    "cores-per-worker": 8,
    "mem-per-worker": "32gb",          # whole job (all cores)
    "walltime": None,                  # e.g. "12:00:00"; None = no request
    "merge-mem": "32gb",
    "scratch-dir": "/tmp",             # job-local working directory (must be bound in the image)
    "apptainer-image": "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-latest",
    "apptainer-binds": ["/cvmfs", "/cms/data", "/home/$USER", "/tmp"],
}


# ---------------------------------------------------------------- helpers
def jobs_dir_of(path):
    """Accept either the pocket-coffea output dir or its job/ subdirectory."""
    for d in (path, os.path.join(path, "job")):
        if os.path.isfile(os.path.join(d, SETTINGS_FILE)):
            return os.path.abspath(d)
    sys.exit(f"No {SETTINGS_FILE} found in {path} or {path}/job")


def load_settings(jobs_dir):
    with open(os.path.join(jobs_dir, SETTINGS_FILE)) as f:
        return json.load(f)


def apptainer_prefix(settings):
    binds = " ".join(f"-B {b}" for b in settings["apptainer-binds"])
    return f'apptainer exec {binds} --env KRB5CCNAME="${{KRB5CCNAME:-}}" {settings["apptainer-image"]}'


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


def qsub(script, depend=None):
    if shutil.which("qsub") is None:
        sys.exit("qsub not found. Inside the apptainer image PBS is not available: run\n"
                 f"    python {os.path.abspath(__file__)} submit <outputdir>\n"
                 "on the host (outside the container).")
    cmd = ["qsub"] + (["-W", f"depend={depend}"] if depend else []) + [script]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
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
cp "$JOBS_DIR/config_job_$JOBID.pkl" "$JOBS_DIR/job.sh" "$JOBS_DIR/{settings['inner-yaml']}" .

{apptainer_prefix(settings)} \\
    bash job.sh "$JOBID" "config_job_$JOBID.pkl" "$CHUNKSIZE" > "$JOBS_DIR/logs/job_$JOBID.log" 2>&1
STATUS=$?

cd / && rm -rf "$WORK"
exit $STATUS
"""
    path = os.path.join(jobs_dir, f"pbs_array_{tag}.sh")
    with open(path, "w") as f:
        f.write(script)
    return path


def write_merge_script(jobs_dir, settings, tag):
    mem_gb = _mem_gb(settings["merge-mem"])
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
    for i in ids:
        for state in ("failed", "running", "done"):
            p = os.path.join(jobs_dir, f"job_{i}.{state}")
            if os.path.exists(p):
                os.remove(p)
        open(os.path.join(jobs_dir, f"job_{i}.idle"), "w").close()

    joblist = write_joblist(jobs_dir, f"joblist_{tag}.txt", ids, settings)
    array_id = qsub(write_array_script(jobs_dir, settings, joblist, tag))
    merge_id = qsub(write_merge_script(jobs_dir, settings, tag), depend=f"afterany:{array_id}")
    record_submission(jobs_dir, {"tag": tag, "array": array_id, "merge": merge_id, "jobs": list(ids)})
    print(f"{settings['name']}: submitted {len(ids)} job(s) as {array_id}; merge job {merge_id} "
          f"-> {settings['output-dir']}/output_all.coffea")
    return array_id, merge_id


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
    merged = os.path.exists(os.path.join(settings["output-dir"], "output_all.coffea"))
    print(f"{settings['name']:24} " + "  ".join(f"{k}={v}" for k, v in counts.items())
          + f"  merged={'yes' if merged else 'no'}")
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


def merge(path):
    jobs_dir = jobs_dir_of(path)
    settings = load_settings(jobs_dir)
    merge_id = qsub(write_merge_script(jobs_dir, settings, "manual"))
    record_submission(jobs_dir, {"tag": "manual-merge", "merge": merge_id})
    print(f"{settings['name']}: merge job {merge_id}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["submit", "status", "resubmit", "merge"])
    ap.add_argument("outputdirs", nargs="+", help="pocket-coffea output dir(s) (or their job/ subdir)")
    ap.add_argument("--include-unfinished", action="store_true",
                    help="resubmit: also jobs left running/idle/missing (e.g. killed by the batch system)")
    args = ap.parse_args()
    for d in args.outputdirs:
        if args.command == "submit":
            submit(d)
        elif args.command == "status":
            status(d)
        elif args.command == "resubmit":
            resubmit(d, args.include_unfinished)
        else:
            merge(d)
