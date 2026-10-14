#!/usr/bin/env python3
"""
Driver for the Run 3 (2024) processing:

    inference    pocket-coffea run --cfg inference_conf_<group>.py   (per group, era set, jet variation)
    status       PBS backend: state of the jobs and merges
    resubmit     PBS backend: resubmit failed jobs (+ a new merge)

Pickling and datacards for Run 3 are not ported yet (steps after the trigger, b-tag and Xbb SFs).

Outputs are organised as
    <outdir>/nominal/<group>/output_all.coffea
    <outdir>/<JEC variation>/<group>/output_all.coffea

Examples:
    python run_pipeline.py inference --backend pbs --groups trigeff
    python run_pipeline.py inference --backend pbs --groups signal bkg minor data
    python run_pipeline.py status --backend pbs --groups signal bkg minor data
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
STEPS = ['inference', 'status', 'resubmit']
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
    ap.add_argument('--eras', nargs='+', default=['2024'])
    ap.add_argument('--jec-variations', nargs='*', default=[],
                    help='jet variations to run in addition to nominal (PC_JEC_VARIATION values)')
    ap.add_argument('--outdir', default='outputs')
    ap.add_argument('--backend', choices=['local', 'pbs'], default='local')
    ap.add_argument('--container', choices=['auto', 'always', 'never'], default='auto',
                    help='wrap commands in apptainer exec (auto: when not already in a container)')
    ap.add_argument('--include-unfinished', action='store_true',
                    help='resubmit: also jobs left running/idle (e.g. killed by the batch system)')
    ap.add_argument('--dry-run', action='store_true', help='print the commands only')
    args = ap.parse_args(argv)
    args.pc_args = pc_args

    steps = STEPS if 'all' in args.steps else [s for s in STEPS if s in args.steps]
    if args.backend == 'pbs' and 'inference' in steps and len(steps) > 1:
        print("PBS backend: the inference jobs run asynchronously; only submitting them now.\n"
              "Run `status` / `resubmit` afterwards.")
        steps = ['inference']
    for step in steps:
        print(f"\n===================== {step} =====================")
        globals()[f"step_{step}"](args)
