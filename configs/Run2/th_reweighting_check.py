#!/usr/bin/env python3
"""
How much does the tH SM reweighting matter? Run 2 check for the Run 3 decision (the 2024 THQ / THW
ctcvcp samples are generated at kappa_t = -1 like the UL ones, but their NanoAOD has no
LHEReweightingWeight, so in Run 3 they can only be used as they are).

Per tH dataset (THQ_*, THW_*) of the .coffea outputs, two yields with the same SM cross section and lumi:

    SM reweighted (as the Run 2 cards):  w * th_sm_weight * L sigma / sum_all sign(gw) * th_sm_weight  (theory_sums)
    as is (kappa_t = -1 kinematics):     w                * L sigma / sum_all sign(gw)                (sum_signOf_genweights)

w = sign(genWeight) x the other nominal weight factors (outvars.nominal_weight_factors without th_sm_weight).
Both have the right total rate; they differ only in acceptance and shapes. Events: makeDatacard.cuts,
disc_<model> >= 0, binned as the default cards (reco pT channels Zhpt1-3, NN edges of nn_edges_disc_<model>.json,
summed over the mass bins). The total background per bin (TTBar, tt_B, QCD, VJets, SingleTop, ttX) comes from
the nominal pickles, so the difference can be compared with it.

    python3 th_reweighting_check.py -i outputs/nominal/minor [--btag-norm btag_norm.json] [--model M]

-i takes .coffea files or output dirs (merged/<dataset>.coffea or output_all.coffea). Pass --btag-norm with
the same JSON as outputToPkl, so that the tH weights match the pickles (the factor is almost common to both).
THQ_2017 BC0B1DEE (68 points, th_sm_weight = 0): left out of the reweighted yield and its normalization sum,
kept in the as-is yield; each is a correct estimate on its own.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from analysis_years import ERA_TO_YEAR, LUMI_FB, PICKLE_DIR, find_era
from makeDatacard import cuts, pt_bins
from outputToPkl import CATEGORY, XSEC_DICT, ERA_SUFFIX
from outvars import nominal_weight_factors
from spanet_models import ACTIVE_MODELS, score_columns
from th_reweighting import is_th

BKG_PROCESSES = ['TTBar', 'tt_B', 'QCD', 'VJets', 'SingleTop', 'ttX']
CHANNEL_EDGES = pt_bins + [500]           # as DataCardShapes: (lo, hi], ZH_pt clipped into [1, 499]
CHANNELS = range(1, len(CHANNEL_EDGES) - 1)
SELECTION_VARS = ['n_ak4jets', 'n_b_outZH', 'ZH_bbvLscore', 'MET_pt', 'ZH_M', 'ZH_pt']


def coffea_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            merged = sorted(glob.glob(os.path.join(p, "merged", "*.coffea")))
            files += merged or sorted(glob.glob(os.path.join(p, "output_all.coffea")))
        else:
            files.append(p)
    return files


def add_disc(df, model):
    s, b, l = (df[f"{k}_score_{model}"].to_numpy(dtype=float) for k in ("sig", "ttbb", "ttlf"))
    den = s + b + l
    df["disc"] = np.where(den > 0, s / np.where(den > 0, den, 1.0), 0.0)
    return df


def binned(df, weight, edges, year):
    """{channel: (sum w, sum w^2) per NN bin} after the card selection."""
    sel = df[cuts(df, year) & (df["disc"] >= 0)]
    pt = sel["ZH_pt"].clip(CHANNEL_EDGES[0] + 1, CHANNEL_EDGES[-1] - 1)
    ch = pd.cut(pt, bins=CHANNEL_EDGES, labels=False)
    out = {}
    for i in CHANNELS:
        m = (ch == i).to_numpy()
        x, w = sel["disc"].to_numpy()[m], sel[weight].to_numpy()[m]
        e = np.asarray(edges[f"Zhpt{i}"], dtype=float)
        out[i] = (np.histogram(x, bins=e, weights=w)[0], np.histogram(x, bins=e, weights=w * w)[0])
    return out


def add_to(acc, key, hists):
    for i, (w, w2) in hists.items():
        if (key, i) in acc:
            acc[(key, i)] = (acc[(key, i)][0] + w, acc[(key, i)][1] + w2)
        else:
            acc[(key, i)] = (w.copy(), w2.copy())


def read_th(files, model, btag_norm):
    variables = ['genWeight', 'th_sm_weight'] + SELECTION_VARS + score_columns(model) + list(nominal_weight_factors)
    frames, norms = [], {}
    for path in files:
        from coffea.util import load
        out = load(path)
        for sample, by_ds in out["columns"].items():
            for dataset, by_cat in by_ds.items():
                if not is_th(dataset) or "nominal" not in by_cat.get(CATEGORY, {}):
                    continue
                cols = by_cat[CATEGORY]["nominal"]
                df = pd.DataFrame({v: cols[f"events_{v}"].value for v in dict.fromkeys(variables)
                                   if f"events_{v}" in cols})
                era = find_era(dataset)
                if btag_norm is not None:
                    from btag_norm import apply_btag_norm
                    df = apply_btag_norm(df, era, sample, btag_norm, "nominal")
                base = sample.split("__")[0]
                lumi_xs = LUMI_FB[era] * 1000.0 * XSEC_DICT[ERA_SUFFIX.sub("", dataset)]
                norms[dataset] = (lumi_xs / out["theory_sums"][base][dataset]["all"]["nominal"],
                                  lumi_xs / out["sum_signOf_genweights"][dataset])
                w = np.sign(df["genWeight"].to_numpy())
                for c in nominal_weight_factors:
                    if c != "th_sm_weight" and c in df:
                        w = w * df[c].to_numpy()
                df["w_rw"] = w * df["th_sm_weight"].to_numpy() * norms[dataset][0]
                df["w_asis"] = w * norms[dataset][1]
                df["dataset"], df["era"], df["year"] = dataset, era, ERA_TO_YEAR[era]
                frames.append(add_disc(df, model))
                print(f"  {dataset}: {len(df)} events, norm SM-rw {norms[dataset][0]:.4g}, as is {norms[dataset][1]:.4g}")
        del out
    if not frames:
        raise SystemExit("No tH datasets (THQ_*, THW_*) in the inputs")
    return pd.concat(frames, ignore_index=True)


def background(year, model, edges):
    acc = {}
    for proc in BKG_PROCESSES:
        path = os.path.join(PICKLE_DIR.format(year=year), f"{proc}.pkl")
        if not os.path.exists(path):
            print(f"  [WARNING] {path} missing: left out of the background")
            continue
        from datacard_shapes import getZhbbWeight
        df = pd.read_pickle(path)
        df["disc"] = df[f"disc_{model}"] if f"disc_{model}" in df else add_disc(df, model)["disc"]
        df["_w"] = getZhbbWeight(df, year)
        add_to(acc, "bkg", binned(df, "_w", edges, year))
    return {i: acc[("bkg", i)][0] for i in CHANNELS} if acc else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True, help="minor-group .coffea files or output dirs")
    ap.add_argument("--model", default=ACTIVE_MODELS[0])
    ap.add_argument("--edges", default=None, help="default nn_edges_disc_<model>.json")
    ap.add_argument("--btag-norm", default=None, help="btag_norm.json (as outputToPkl)")
    args = ap.parse_args()

    import json
    edges_all = json.load(open(args.edges or f"nn_edges_disc_{args.model}.json"))
    btag_norm = json.load(open(args.btag_norm)) if args.btag_norm else None
    th = read_th(coffea_files(args.input), args.model, btag_norm)

    print(f"\nModel {args.model}; selection makeDatacard.cuts, disc >= 0; NN bins of the default cards\n")
    print("Yield after the selection (SM-reweighted vs as is, same SM cross section):")
    for (year, proc), d in th.assign(proc=th["dataset"].str.split("_").str[0]).groupby(["year", "proc"]):
        sel = d[cuts(d, year) & (d["disc"] >= 0)]
        rw, asis = sel["w_rw"].sum(), sel["w_asis"].sum()
        print(f"  {year} {proc}: SM-rw {rw:8.3f} +- {np.sqrt((sel['w_rw'] ** 2).sum()):.3f}   as is {asis:8.3f} "
              f"+- {np.sqrt((sel['w_asis'] ** 2).sum()):.3f}   as is / SM-rw {asis / rw if rw else float('nan'):.3f}")

    for year, d in th.groupby("year"):
        edges = edges_all.get(year, next(iter(edges_all.values())))
        acc = {}
        add_to(acc, "rw", binned(d, "w_rw", edges, year))
        add_to(acc, "asis", binned(d, "w_asis", edges, year))
        bkg = background(year, args.model, edges)
        print(f"\n{year}: tH (THQ + THW) per NN bin; diff = as is - SM-rw; % of the total background")
        print(f"  {'bin':9s} {'SM-rw':>9s} {'as is':>9s} {'ratio':>6s} {'bkg':>10s} {'tH/bkg':>7s} {'diff/bkg':>9s}")
        worst = 0.0
        for i in CHANNELS:
            (rw, rw2), (asis, _) = acc[("rw", i)], acc[("asis", i)]
            for k in range(len(rw)):
                b = bkg[i][k] if bkg is not None else float("nan")
                share, diff = rw[k] / b if b else float("nan"), (asis[k] - rw[k]) / b if b else float("nan")
                worst = max(worst, abs(diff)) if np.isfinite(diff) else worst
                ratio = asis[k] / rw[k] if rw[k] else float("nan")
                print(f"  Zhpt{i} NN{k + 1:<2d} {rw[k]:9.4f} {asis[k]:9.4f} {ratio:6.2f} {b:10.2f} "
                      f"{100 * share:6.3f}% {100 * diff:+8.3f}%   (SM-rw stat {np.sqrt(rw2[k]) / rw[k] if rw[k] else 0:.0%})")
        print(f"  largest |as is - SM-rw| / bkg: {100 * worst:.3f}%")


if __name__ == "__main__":
    main()
