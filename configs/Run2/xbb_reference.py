#!/usr/bin/env python3
"""
Run 2 reference for the Run 3 Xbb working point: process composition after the full analysis
selection (makeDatacard.cuts, including Xbb >= LP of the era, 0.9172 in 2018), from the pickles.

The Run 2 outputs only hold events above the LP (the processor applies it), so the efficiency of the
Run 2 cut cannot be measured; the composition can. Yields (datacard weight, getZhbbWeight) per class
and pT channel (ZH_pt [200, 300, 450, inf), as Run3/scan_xbb_wp.py):

    ttZ, ttH   signal with gen Z/H pT >= 200 GeV (the measured part)
    sig_lowpt  ttZ / ttH with gen pT < 200 GeV
    tt_B       TTbb samples, tt+B
    tt+C       TTBar with tt_C (genTtbarId % 100 in 41..49)
    tt+LF      the rest of TTBar
    minor      SingleTop, ttX, VJets, QCD

Prints the yields and the ratios to the signal (ttZ + ttH, and ttH alone); writes a JSON for
Run3/scan_xbb_wp.py composition --ref, which finds the Run 3 cut with the same ratios.

    python3 xbb_reference.py --year 2018                  # -> xbb_reference_2018.json

Caveats: the ttZ normalization is the cross section the pickles were made with (+30% since 4222ca8;
use the ttH ratios if the pickles are older); the Run 2 weight includes the Xbb SF (xbb_sf), Run 3 has
none yet; the b tagger (DeepJet vs UParT) and sqrt(s) also change the composition, not only the Xbb cut.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from analysis_years import PICKLE_DIR, bbvl_wp
from datacard_shapes import getZhbbWeight
from makeDatacard import cuts
from outvars import nominal_weight_factors

PT_EDGES = [200, 300, 450, np.inf]
SIG_GENPT = 200.0
CLASSES = ["ttZ", "ttH", "sig_lowpt", "tt_B", "tt+C", "tt+LF", "minor"]
PICKLES = {"ttZ": "ttZ", "ttH": "ttH", "tt_B": "tt_B", "TTBar": "TTBar",
           "SingleTop": "minor", "ttX": "minor", "VJets": "minor", "QCD": "minor"}
READ = ["ZH_pt", "ZH_M", "MET_pt", "n_ak4jets", "n_b_outZH", "ZH_bbvLscore", "era", "norm_weight",
        "genWeight", "genZHpt", "tt_C"] + list(nominal_weight_factors)


def class_of(pickle, df):
    """Class label per event of one pickle."""
    lab = PICKLES[pickle]
    if lab in ("ttZ", "ttH"):
        return np.where(df["genZHpt"].to_numpy() >= SIG_GENPT, lab, "sig_lowpt")
    if lab == "TTBar":
        return np.where(df["tt_C"].to_numpy(dtype=bool), "tt+C", "tt+LF")
    return np.full(len(df), lab)


def composition(year):
    """{'total': {class: yield}, 'channels': [{class: yield} per pT channel]}"""
    tot = {c: 0.0 for c in CLASSES}
    chans = [{c: 0.0 for c in CLASSES} for _ in PT_EDGES[:-1]]
    for pickle in PICKLES:
        path = f"{PICKLE_DIR.format(year=year)}/{pickle}.pkl"
        if not os.path.exists(path):
            print(f"  [warning] {path} not found: skipped")
            continue
        df = pd.read_pickle(path)
        df = df[[c for c in READ if c in df]]
        df = df[cuts(df, year) & (df["ZH_pt"] >= PT_EDGES[0])]
        w = getZhbbWeight(df, year).to_numpy(dtype=float)
        cls = class_of(pickle, df)
        ch = np.digitize(df["ZH_pt"].to_numpy(), PT_EDGES) - 1
        for c in np.unique(cls):
            m = cls == c
            tot[c] += float(w[m].sum())
            for i in range(len(chans)):
                chans[i][c] += float(w[m & (ch == i)].sum())
        print(f"  {pickle}: {len(df)} events after the selection")
    return {"total": tot, "channels": chans}


def print_table(comp, title):
    """Yields and ratios to the signal, total and per pT channel (shared with Run3/scan_xbb_wp.py)."""
    rows = [("total", comp["total"])] + [
        (f"{PT_EDGES[i]:.0f}-{PT_EDGES[i + 1]:.0f}", c) for i, c in enumerate(comp["channels"])]
    print(f"\n{title}")
    print(f"{'pT':>9} " + " ".join(f"{c:>10}" for c in CLASSES)
          + f" {'LF/S':>7} {'C/S':>7} {'B/S':>7} {'LF/ttH':>7} {'C/ttH':>7} {'B/ttH':>7}")
    for name, y in rows:
        s, h = y["ttZ"] + y["ttH"], y["ttH"]
        r = lambda a, b: a / b if b > 0 else np.nan
        print(f"{name:>9} " + " ".join(f"{y[c]:10.1f}" for c in CLASSES)
              + f" {r(y['tt+LF'], s):7.2f} {r(y['tt+C'], s):7.2f} {r(y['tt_B'], s):7.2f}"
              + f" {r(y['tt+LF'], h):7.2f} {r(y['tt+C'], h):7.2f} {r(y['tt_B'], h):7.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", default="2018", help="pickle year (PICKLE_DIR), e.g. 2018, 2017, 2016")
    ap.add_argument("-o", "--output", help="default: xbb_reference_<year>.json")
    args = ap.parse_args()
    comp = composition(args.year)
    era = args.year if args.year != "2016" else "2016_PostVFP"
    print_table(comp, f"Run 2 {args.year}: composition after the full selection (Xbb >= LP, {bbvl_wp(era)} in {era})")
    out = args.output or f"xbb_reference_{args.year}.json"
    with open(out, "w") as f:
        json.dump({"year": args.year, "pt_edges": [e if np.isfinite(e) else None for e in PT_EDGES],
                   "classes": CLASSES, **comp}, f, indent=1)
    print(f"\nWrote {out} (for Run3/scan_xbb_wp.py composition --ref {os.path.abspath(out)})")
