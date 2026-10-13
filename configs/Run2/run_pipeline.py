#!/usr/bin/env python3
"""
Driver for the analysis chain after the SPANet training:

    inference    pocket-coffea run --cfg inference_conf_<group>.py   (per group, era set, jet variation)
    btagnorm     btag_norm.py on the nominal outputs                  -> btag_norm.json
    pickle       outputToPkl.py (nominal + jet-variation outputs)     -> pickled/Nominal<year>/
    datacard     makeDatacard.py

Outputs are organised as
    <outdir>/nominal/<group>/output_all.coffea
    <outdir>/<JEC variation>/<group>/output_all.coffea

Examples:
    python run_pipeline.py inference --groups signal bkg minor --eras 2018
    python run_pipeline.py inference --groups signal bkg minor --jec-variations JES_Total JER
    python run_pipeline.py btagnorm pickle datacard
    python run_pipeline.py all --dry-run          # print every command without running it

Anything after '--' is passed on to 'pocket-coffea run', e.g.  -- --executor futures --test
"""
import argparse
import os
import shlex
import subprocess
import sys

from samples import GROUPS

HERE = os.path.dirname(os.path.abspath(__file__))
MC_GROUPS = ['signal', 'bkg', 'minor']
STEPS = ['inference', 'btagnorm', 'pickle', 'datacard']


def run(cmd, env_extra=None, dry=False):
    env_str = " ".join(f"{k}={v}" for k, v in (env_extra or {}).items())
    print(f"\n>>> {env_str + ' ' if env_str else ''}{' '.join(shlex.quote(c) for c in cmd)}", flush=True)
    if dry:
        return
    env = dict(os.environ, **(env_extra or {}))
    subprocess.run(cmd, cwd=HERE, env=env, check=True)


def coffea_file(outdir, tag, group):
    return os.path.join(outdir, tag, group, "output_all.coffea")


def step_inference(args):
    for tag in ['nominal'] + args.jec_variations:
        for group in args.groups:
            if tag != 'nominal' and GROUPS[group].get('is_data'):
                continue                      # no jet variations for data
            env = {"PC_ERAS": ",".join(args.eras), "PC_JEC_VARIATION": tag}
            cmd = ["pocket-coffea", "run", "--cfg", f"inference_conf_{group}.py",
                   "-o", os.path.join(args.outdir, tag, group)] + args.pc_args
            run(cmd, env, args.dry_run)


def existing(paths, dry):
    found = [p for p in paths if dry or os.path.exists(p)]
    for p in set(paths) - set(found):
        print(f"  [skip] {p} not found")
    return found


def step_btagnorm(args):
    # nominal + jet-variation outputs: r is derived per shape variation (JES-varied b-tag SF)
    inputs = existing([coffea_file(args.outdir, tag, g) for tag in ['nominal'] + args.jec_variations
                       for g in args.groups if not GROUPS[g].get('is_data')], args.dry_run)
    if inputs:
        run([sys.executable, "btag_norm.py", "-i", *inputs, "-o", args.btag_norm], dry=args.dry_run)


def step_pickle(args):
    mc_groups = [g for g in args.groups if not GROUPS[g].get('is_data')]
    nominal = existing([coffea_file(args.outdir, 'nominal', g) for g in mc_groups], args.dry_run)
    btag = ["--btag-norm", args.btag_norm] if (args.dry_run or os.path.exists(args.btag_norm)) else []
    if nominal:
        run([sys.executable, "outputToPkl.py", "-i", *nominal, *btag], dry=args.dry_run)
    for tag in args.jec_variations:
        varied = existing([coffea_file(args.outdir, tag, g) for g in mc_groups], args.dry_run)
        if varied:
            run([sys.executable, "outputToPkl.py", "-i", *varied, *btag,
                 "--variations", "all", "--skip-nominal"], dry=args.dry_run)


def step_datacard(args):
    run([sys.executable, "makeDatacard.py"] + (["--eft-ttbb"] if args.eft_ttbb else []), dry=args.dry_run)


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
    ap.add_argument('--btag-norm', default='btag_norm.json')
    ap.add_argument('--eft-ttbb', action='store_true', help='passed to makeDatacard.py')
    ap.add_argument('--dry-run', action='store_true', help='print the commands only')
    args = ap.parse_args(argv)
    args.pc_args = pc_args

    steps = STEPS if 'all' in args.steps else [s for s in STEPS if s in args.steps]
    for step in steps:
        print(f"\n===================== {step} =====================")
        globals()[f"step_{step}"](args)
