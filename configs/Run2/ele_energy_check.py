#!/usr/bin/env python3
"""
Electron energy scale / smearing check over the full analysis selection (Run 2 UL, NanoAODv9).

EGM (UL): the NanoAODv9 electron pT already carries the scale (data) and smearing (MC) corrections;
the uncertainties exist as Electron_dEscaleUp/Down (data only) and Electron_dEsigmaUp/Down (MC only).
The inference stores them for the leading selected electron as relative shifts dE/E
(outvars.ele_energy_vars, get_controlvars). This script shifts the electron pT by (1 + dE/E) and
reports the fraction of the selected yield (makeDatacard.cuts, datacard weight) whose electron
falls below the selection threshold (ELE_PT_MIN), per process and Z/H pT channel:

    MC   (smearing, dEsigma): from the pickles (PICKLE_DIR / <process>.pkl)
    data (scale, dEscale):    from the data .coffea output(s) with --data (unweighted counts)

Only migration out of the selection is visible (events just below the threshold are not stored);
migration in is about the mirror image of the opposite variation.

    python3 ele_energy_check.py
    python3 ele_energy_check.py --years 2018 --data output_data.coffea
"""
import argparse
import os

import numpy as np
import pandas as pd

from analysis_years import YEARS, PICKLE_DIR, find_era, ERA_TO_YEAR
from datacard_shapes import getZhbbWeight

ELE_PT_MIN = 35.0          # params/skim_params.yaml Electron pt (lep_sel: pt > cut)
CHANNELS = {"Zhpt1": (200, 300), "Zhpt2": (300, 450), "Zhpt3": (450, np.inf)}
PROCESSES = ["ttZ", "ttH", "tt_B", "TTBar", "SingleTop", "ttX", "VJets"]
SHIFTS = {"MC": ["ele_dEsigmaUp", "ele_dEsigmaDown"], "data": ["ele_dEscaleUp", "ele_dEscaleDown"]}
NEEDED = ["ele_pt", "ZH_pt"]


def migration(df, w, shift):
    """Fraction of the selected yield lost (all events, electron events) for one shift column."""
    pt = df["ele_pt"].to_numpy(dtype=float)
    ele = np.nan_to_num(pt, nan=0.0) > 0
    rel = np.nan_to_num(df[shift].to_numpy(dtype=float), nan=0.0)
    lost = ele & (pt * (1 + rel) <= ELE_PT_MIN)
    tot, tot_e = w.sum(), w[ele].sum()
    return (w[lost].sum() / tot if tot else 0.0), (w[lost].sum() / tot_e if tot_e else 0.0)


def report(label, df, w, shifts):
    rows = []
    for ch, (lo, hi) in {"all": (-np.inf, np.inf), **CHANNELS}.items():
        m = ((df["ZH_pt"] > lo) & (df["ZH_pt"] <= hi)).to_numpy()
        if not m.any():
            continue
        for s in shifts:
            f_all, f_ele = migration(df[m], w[m], s)
            rel = np.abs(np.nan_to_num(df.loc[m, s].to_numpy(dtype=float)))
            rows.append((label, ch, s.replace("ele_", ""), f_all, f_ele, np.percentile(rel[rel > 0], 99) if (rel > 0).any() else 0.0))
    for r in rows:
        print(f"{r[0]:>12} {r[1]:>6} {r[2]:>12}  lost {r[3]:9.2e} of all, {r[4]:9.2e} of e events;  99% |dE/E| {r[5]:.1e}")


def check_mc(years):
    from makeDatacard import cuts
    for year in years:
        print(f"--- MC {year} (smearing, datacard weight)")
        for proc in PROCESSES:
            path = os.path.join(PICKLE_DIR.format(year=year), f"{proc}.pkl")
            if not os.path.exists(path):
                continue
            df = pd.read_pickle(path)
            missing = [c for c in NEEDED + SHIFTS["MC"] if c not in df]
            if missing:
                print(f"{proc:>12}: no {missing} (pickles made before the ele_energy_vars columns)")
                continue
            df = df[cuts(df, year)]
            report(proc, df, getZhbbWeight(df, year).to_numpy(dtype=float), SHIFTS["MC"])


def data_frame(paths, years):
    """Selected-category columns of the data .coffea outputs as one DataFrame per year."""
    from coffea.util import load
    cols = NEEDED + SHIFTS["data"] + ["n_ak4jets", "n_b_outZH", "ZH_bbvLscore", "MET_pt", "ZH_M"]
    frames = {}
    for path in paths:
        out = load(path)
        for sample, datasets in out["columns"].items():
            for dataset, cats in datasets.items():
                c = cats.get("btag_mask", {}).get("nominal")
                if c is None:
                    continue
                era = find_era(dataset)
                year = ERA_TO_YEAR[era]
                if year not in years:
                    continue
                if any(f"events_{v}" not in c for v in cols):
                    print(f"  {dataset}: missing columns, skipped")
                    continue
                df = pd.DataFrame({v: c[f"events_{v}"].value for v in cols})
                df["era"] = era
                frames.setdefault(year, []).append(df)
    return {y: pd.concat(f, ignore_index=True) for y, f in frames.items()}


def check_data(paths, years):
    from makeDatacard import cuts
    for year, df in sorted(data_frame(paths, years).items()):
        print(f"--- data {year} (scale, unweighted)")
        df = df[cuts(df, year)]
        report("data", df, np.ones(len(df)), SHIFTS["data"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", nargs="+", default=list(YEARS))
    ap.add_argument("--data", nargs="*", default=[], help="data .coffea output(s) for the scale check")
    ap.add_argument("--no-mc", action="store_true")
    args = ap.parse_args()
    if not args.no_mc:
        check_mc(args.years)
    if args.data:
        check_data(args.data, args.years)


if __name__ == "__main__":
    main()
