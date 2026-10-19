#!/usr/bin/env python3
"""
Driver for the analysis chain after the SPANet training:

    inference    pocket-coffea run --cfg inference_conf_<group>.py   (per group, era set, jet variation)
    status       PBS backend: state of the jobs and merges
    resubmit     PBS backend: resubmit failed jobs (+ a new merge)
    btagnorm     btag_norm.py on the nominal + variation outputs     -> btag_norm.json
    pickle       outputToPkl.py (nominal + jet-variation outputs)     -> pickled/Nominal<year>/
    norms        make_process_norms.py (nominal outputs)              -> process_norms/process_norms_<year>.json
    datacard     makeDatacard.py

Outputs are organised as
    <outdir>/nominal/<group>/output_all.coffea
    <outdir>/<JEC variation>/<group>/output_all.coffea

Backends for the inference step:
    local  pocket-coffea runs here (executor from the config / arguments after '--')
    pbs    one PBS Pro job array per group and variation (../common/executors_pbs.py, run_options_pbs.yaml):
           jobs are prepared in the apptainer image, submitted from the host, and merged by a
           PBS job into output_all.coffea. Follow with `status` / `resubmit`; run btagnorm,
           pickle and datacard once everything is merged.

Run this script on the host or inside the apptainer image: with --container auto (default),
pocket-coffea and the python steps are wrapped in `apptainer exec` when not already in a container.

Examples:
    python run_pipeline.py inference --groups signal bkg minor --eras 2018 --backend pbs
    python run_pipeline.py inference --backend pbs --jec-variations JES_Total JER
    python run_pipeline.py status --backend pbs --jec-variations JES_Total JER
    python run_pipeline.py resubmit --backend pbs
    python run_pipeline.py btagnorm pickle datacard --jec-variations JES_Total JER
    python run_pipeline.py inference --dry-run      # print every command without running it

Anything after '--' is passed on to 'pocket-coffea run', e.g.  -- --executor futures --test
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
COMMON = os.path.join(os.path.dirname(HERE), "common")       # shared PBS tools (configs/common)
sys.path.insert(0, COMMON)

import pbs_tools  # noqa: E402
from samples import GROUPS  # noqa: E402

MC_GROUPS = ['signal', 'bkg', 'minor']
STEPS = ['inference', 'status', 'resubmit', 'btagnorm', 'pickle', 'norms', 'datacard']
PBS_RUN_OPTIONS = "run_options_pbs.yaml"


# ---------------------------------------------------------------- execution helpers
def in_container():
    return bool(os.environ.get("APPTAINER_CONTAINER") or os.environ.get("SINGULARITY_CONTAINER"))


def container_settings():
    """Image and binds: run_options_pbs.yaml if PyYAML is available, else pbs_tools defaults."""
    settings = {k: pbs_tools.DEFAULTS[k] for k in ("apptainer-image", "apptainer-binds")}
    try:
        import yaml
        with open(os.path.join(HERE, PBS_RUN_OPTIONS)) as f:
            opts = yaml.safe_load(f) or {}
        settings.update({k: opts[k] for k in settings if opts.get(k)})
    except ImportError:
        pass
    return settings


def wrap(cmd, args):
    use = args.container == "always" or (args.container == "auto" and not in_container())
    if not use:
        return cmd
    s = container_settings()
    binds = [x for b in s["apptainer-binds"] for x in ("-B", os.path.expandvars(b))]
    return (["apptainer", "exec"] + binds + ["--env", f"KRB5CCNAME={os.environ.get('KRB5CCNAME', '')}",
                                             s["apptainer-image"]] + cmd)


def run(cmd, args, env_extra=None):
    cmd = wrap(cmd, args)
    env_str = " ".join(f"{k}={v}" for k, v in (env_extra or {}).items())
    print(f"\n>>> {env_str + ' ' if env_str else ''}{' '.join(shlex.quote(c) for c in cmd)}", flush=True)
    if args.dry_run:
        return
    env = dict(os.environ, **(env_extra or {}))
    subprocess.run(cmd, cwd=HERE, env=env, check=True)


def python(script, *script_args):
    return ["python3", script, *script_args]


def runs(args, include_data=True):
    """(tag, group, output dir) of every inference run."""
    for tag in ['nominal'] + args.jec_variations:
        for group in args.groups:
            if GROUPS[group].get('is_data') and (tag != 'nominal' or not include_data):
                continue                      # no jet variations for data
            yield tag, group, os.path.join(args.outdir, tag, group)


def coffea_file(outdir, tag, group):
    return os.path.join(outdir, tag, group, "output_all.coffea")


def existing(paths, dry):
    found = [p for p in paths if dry or os.path.exists(p)]
    for p in sorted(set(paths) - set(found)):
        print(f"  [skip] {p} not found")
    return found


# ---------------------------------------------------------------- steps
def step_inference(args):
    can_submit = shutil.which("qsub") is not None
    for tag, group, outdir in runs(args):
        env = {"PC_ERAS": ",".join(args.eras), "PC_JEC_VARIATION": tag}
        cmd = ["pocket-coffea", "run", "--cfg", f"inference_conf_{group}.py", "-o", outdir]
        if args.backend == "pbs":
            cmd += ["--executor", "pbs", "--executor-custom-setup", os.path.join(COMMON, "executors_pbs.py"),
                    "-ro", PBS_RUN_OPTIONS]
            env[pbs_tools.NO_SUBMIT_ENV] = "1"       # the jobs are submitted below, from the host
        run(cmd + args.pc_args, args, env)

        if args.backend == "pbs":
            if args.dry_run:
                print(f">>> pbs_tools.submit({outdir})")
            elif can_submit:
                pbs_tools.submit(os.path.join(HERE, outdir))
            else:
                print(f"  qsub not available here: submit from the host with\n"
                      f"    python {pbs_tools.__file__} submit {os.path.join(HERE, outdir)}")


def _pbs_dirs(args):
    for tag, group, outdir in runs(args):
        path = os.path.join(HERE, outdir)
        if os.path.isfile(os.path.join(path, "job", pbs_tools.SETTINGS_FILE)):
            yield path
        else:
            print(f"  [skip] no PBS jobs in {outdir}")


def step_status(args):
    if args.backend != "pbs":
        sys.exit("status needs --backend pbs")
    for path in _pbs_dirs(args):
        pbs_tools.status(path)


def step_resubmit(args):
    if args.backend != "pbs":
        sys.exit("resubmit needs --backend pbs")
    for path in _pbs_dirs(args):
        if args.dry_run:
            print(f">>> pbs_tools.resubmit({path})")
        else:
            pbs_tools.resubmit(path, include_unfinished=args.include_unfinished)


def step_btagnorm(args):
    # nominal + jet-variation outputs: r is derived per shape variation (JES-varied b-tag SF)
    inputs = existing([coffea_file(args.outdir, tag, g) for tag, g, _ in runs(args, include_data=False)],
                      args.dry_run)
    if inputs:
        run(python("btag_norm.py", "-i", *inputs, "-o", args.btag_norm), args)


def step_pickle(args):
    mc_groups = [g for g in args.groups if not GROUPS[g].get('is_data')]
    nominal = existing([coffea_file(args.outdir, 'nominal', g) for g in mc_groups], args.dry_run)
    btag = ["--btag-norm", args.btag_norm] if (args.dry_run or os.path.exists(args.btag_norm)) else []
    if nominal:
        run(python("outputToPkl.py", "-i", *nominal, *btag), args)
    for tag in args.jec_variations:
        varied = existing([coffea_file(args.outdir, tag, g) for g in mc_groups], args.dry_run)
        if varied:
            run(python("outputToPkl.py", "-i", *varied, *btag, "--variations", "all", "--skip-nominal"), args)


def step_norms(args):
    # factors of the theory variations, from the full-sample sums in the nominal outputs
    mc_groups = [g for g in args.groups if not GROUPS[g].get('is_data')]
    nominal = existing([coffea_file(args.outdir, 'nominal', g) for g in mc_groups], args.dry_run)
    if nominal:
        run(python("make_process_norms.py", "-i", *nominal), args)


def step_datacard(args):
    run(python("makeDatacard.py", *(["--eft-ttbb"] if args.eft_ttbb else [])), args)


if __name__ == '__main__':
    argv = sys.argv[1:]
    pc_args = []
    if '--' in argv:
        i = argv.index('--')
        argv, pc_args = argv[:i], argv[i + 1:]

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('steps', nargs='+', choices=STEPS + ['all'])
    ap.add_argument('--groups', nargs='+', default=MC_GROUPS, choices=list(GROUPS),
                    help=f'sample groups (default: {" ".join(MC_GROUPS)})')
    ap.add_argument('--eras', nargs='+', default=['2018'])
    ap.add_argument('--jec-variations', nargs='*', default=[],
                    help='jet variations to run in addition to nominal (PC_JEC_VARIATION values)')
    ap.add_argument('--outdir', default='outputs')
    ap.add_argument('--backend', choices=['local', 'pbs'], default='local')
    ap.add_argument('--container', choices=['auto', 'always', 'never'], default='auto',
                    help='wrap commands in apptainer exec (auto: when not already in a container)')
    ap.add_argument('--include-unfinished', action='store_true',
                    help='resubmit: also jobs left running/idle (e.g. killed by the batch system)')
    ap.add_argument('--btag-norm', default='btag_norm.json')
    ap.add_argument('--eft-ttbb', action='store_true', help='passed to makeDatacard.py')
    ap.add_argument('--dry-run', action='store_true', help='print the commands only')
    args = ap.parse_args(argv)
    args.pc_args = pc_args

    if 'all' in args.steps:
        steps = ['inference', 'btagnorm', 'pickle', 'datacard']
    else:
        steps = [s for s in STEPS if s in args.steps]
    if args.backend == 'pbs' and 'inference' in steps and len(steps) > 1:
        print("PBS backend: the inference jobs run asynchronously; only submitting them now.\n"
              "Run `status` / `resubmit`, then btagnorm pickle datacard once all outputs are merged.")
        steps = ['inference']
    for step in steps:
        print(f"\n===================== {step} =====================")
        globals()[f"step_{step}"](args)
