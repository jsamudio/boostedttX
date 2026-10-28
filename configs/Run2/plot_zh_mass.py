#!/usr/bin/env python3
"""
Z/H-candidate mass distributions in the highest SPANet bins, per ZH-candidate pT channel: stacked
backgrounds, ttZ and ttH (x10) as lines, MC statistical uncertainty hatched (all years summed).

The events pass the datacard selection (makeDatacard.cuts) and fall in a SPANet bin of the datacard
binning (nn_edges_disc_<model>.json, written by makeDatacard.py), in three selections:
    top     - the highest SPANet bin
    second  - the second-highest SPANet bin
    top2    - both together

    python3 plot_zh_mass.py                                  # first active model, YEARS
    python3 plot_zh_mass.py --model refactor_noOverlap_MultiWP --years 2016 2017 2018 -o plots/zh_mass

Writes <outdir>/zh_mass_<model>_<selection>_Zhpt<i>.pdf (and .png).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd

from analysis_years import YEARS, PICKLE_DIR, LUMI_FB, ERA_TO_YEAR
from datacard_shapes import getZhbbWeight
from spanet_models import ACTIVE_MODELS

# reco ZH-candidate pT channels of the datacards (makeDatacard.pt_bins = [0, 200, 300, 450])
PT_CHANNELS = {1: (200, 300), 2: (300, 450), 3: (450, np.inf)}
MASS_BINS = np.arange(50, 205, 5)
SIGNAL_SCALE = 10

GEN_BB_PT_PLOT_BINS = np.array([0, 50, 100, 150, 200, 300, 400, 600, 1000])
# stacked backgrounds, bottom to top (process names of the pickles)
BACKGROUNDS = [
    ("ttX",       ["ttX"],       r"$t\bar{t}t\bar{t}$, $t\bar{t}\gamma$, $t\bar{t}W$", "#b2182b"),
    ("SingleTop", ["SingleTop"], "single top",                                         "#8c564b"),
    ("tt_B",      ["tt_B"],      r"$t\bar{t}+b\bar{b}$",                               "#4daf4a"),
    ("TTBar",     ["TTBar"],     r"$t\bar{t}$+LF, $t\bar{t}+c\bar{c}$",                "#ff9933"),
    ("VJets",     ["VJets"],     "V+jets",                                             "#66ccee"),
]
SIGNALS = [
    ("ttZ", r"$t\bar{t}Z\times$" + str(SIGNAL_SCALE), "#1f3a93", "--"),
    ("ttH", r"$t\bar{t}H\times$" + str(SIGNAL_SCALE), "#e600e6", "--"),
]
SELECTIONS = {"top": [-1], "second": [-2], "top2": [-2, -1]}   # SPANet bin indices (from the top)
SELECTION_LABEL = {"top": "highest SPANet bin", "second": "second-highest SPANet bin",
                   "top2": "two highest SPANet bins"}


def load(years, model):
    """One DataFrame of all processes and years: process, year, channel, score, mass, weight."""
    from makeDatacard import cuts
    score = f"disc_{model}"
    processes = [p for _, procs, _, _ in BACKGROUNDS for p in procs] + [s for s, _, _, _ in SIGNALS]
    frames = []
    for year in years:
        for proc in processes:
            path = os.path.join(PICKLE_DIR.format(year=year), f"{proc}.pkl")
            if not os.path.exists(path):
                print(f"  [skip] {path} not found")
                continue
            df = pd.read_pickle(path)
            if score not in df:
                raise SystemExit(f"{path} has no column {score}: rerun outputToPkl.py with the model active")
            df = df[cuts(df, year) & (df[score] >= 0)]
            frames.append(pd.DataFrame({
                "process": proc, "year": year, "score": df[score].to_numpy(),
                "pt": df["ZH_pt"].to_numpy(), "mass": df["ZH_M"].to_numpy(),
                "weight": getZhbbWeight(df, year).to_numpy()}))
    if not frames:
        raise SystemExit("No pickles found")
    data = pd.concat(frames, ignore_index=True)
    data["channel"] = 0
    for i, (low, high) in PT_CHANNELS.items():
        data.loc[(data["pt"] >= low) & (data["pt"] < high), "channel"] = i
    return data


def in_bins(data, edges, selection):
    """Mask of the events in the selected SPANet bins of their year and channel."""
    mask = np.zeros(len(data), dtype=bool)
    for (year, channel), idx in data.groupby(["year", "channel"]).groups.items():
        if channel == 0:
            continue
        e = np.asarray(edges[year][f"Zhpt{channel}"], dtype=float)
        score = data.loc[idx, "score"].to_numpy()
        sel = np.zeros(len(idx), dtype=bool)
        for b in SELECTIONS[selection]:
            low, high = e[b - 1], e[b]
            sel |= (score >= low) & ((score < high) | ((b == -1) & (score <= high)))
        mask[data.index.get_indexer(idx)] = sel
    return mask


def hist(df):
    sumw, _ = np.histogram(df["mass"], bins=MASS_BINS, weights=df["weight"])
    sumw2, _ = np.histogram(df["mass"], bins=MASS_BINS, weights=df["weight"] ** 2)
    return sumw, sumw2


def score_label(edges, years, channel, selection):
    """SPANet score range of the selected bins, e.g. 'SPANet score > 0.74'; per year if they differ."""
    idx = SELECTIONS[selection]
    def text(e):
        low, high = e[min(idx) - 1], e[max(idx)]
        return f"> {low:.2f}" if high >= 1 else f"{low:.2f} – {high:.2f}"
    per_year = {y: text(np.asarray(edges[y][f"Zhpt{channel}"], dtype=float)) for y in years}
    if len(set(per_year.values())) == 1:
        return f"SPANet score {next(iter(per_year.values()))}"
    return "SPANet score\n" + "\n".join(f"  {y}: {t}" for y, t in per_year.items())

def plot(df, channel, selection, lumi, out, score_text):
    plt.style.use(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(10, 8))
    stack, labels, colors = [], [], []
    total, total_w2 = np.zeros(len(MASS_BINS) - 1), np.zeros(len(MASS_BINS) - 1)
    for _, procs, label, color in BACKGROUNDS:
        sumw, sumw2 = hist(df[df["process"].isin(procs)])
        stack.append(sumw)
        labels.append(label)
        colors.append(color)
        total += sumw
        total_w2 += sumw2
    hep.histplot(stack, MASS_BINS, stack=True, histtype="fill", label=labels, color=colors, ax=ax)
    err = np.sqrt(total_w2)
    ax.stairs(total + err, MASS_BINS, baseline=total - err, fill=True, hatch="xxxx", facecolor="none",
              edgecolor="#d62728", linewidth=0, label="Stat. unc.")
    ymax = np.max(total + err)
    for proc, label, color, ls in SIGNALS:
        sumw, _ = hist(df[df["process"] == proc])
        hep.histplot(SIGNAL_SCALE * sumw, MASS_BINS, histtype="step", color=color, linestyle=ls,
                     linewidth=2, label=label, ax=ax)
        ymax = max(ymax, np.max(SIGNAL_SCALE * sumw))

    low, high = PT_CHANNELS[channel]
    pt_text = (rf"${low} < p_{{T}}^{{Z/H\ cand.}} < {high}$ GeV" if np.isfinite(high)
               else rf"$p_{{T}}^{{Z/H\ cand.}} > {low}$ GeV")
    ax.text(0.04, 0.95, f"{pt_text}\n{score_text}", transform=ax.transAxes,
            ha="left", va="top", fontsize="x-small", linespacing=1.3)
    ax.set_xlim(MASS_BINS[0], MASS_BINS[-1])
    n_lines = 2 + score_text.count("\n")
    ax.set_ylim(0, (1.5 + 0.12 * n_lines) * max(ymax, 1e-9))
    ax.set_xlabel(r"$m_{reg}^{Z/H\ cand.}$ [GeV]")
    ax.set_ylabel(f"Events / {MASS_BINS[1] - MASS_BINS[0]:g} GeV")
    ax.legend(ncol=1, fontsize="x-small", loc="upper right")
    hep.cms.label("Preliminary", data=False, lumi=round(lumi), com=13, ax=ax)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=ACTIVE_MODELS[0] if ACTIVE_MODELS else None)
    ap.add_argument("--years", nargs="+", default=list(YEARS))
    ap.add_argument("--edges", default=None, help="NN edges JSON (default: nn_edges_disc_<model>.json)")
    ap.add_argument("-o", "--outdir", default="plots/zh_mass")
    args = ap.parse_args()

    edges_path = args.edges or f"nn_edges_disc_{args.model}.json"
    if not os.path.exists(edges_path):
        raise SystemExit(f"{edges_path} not found: run makeDatacard.py (it writes the SPANet bin edges)")
    edges = json.load(open(edges_path))
    missing = [y for y in args.years if y not in edges]
    if missing:
        raise SystemExit(f"No bin edges for {missing} in {edges_path}")

    lumi = sum(l for era, l in LUMI_FB.items() if ERA_TO_YEAR[era] in args.years)
    data = load(args.years, args.model)
    for selection in SELECTIONS:
        selected = data[in_bins(data, edges, selection)]
        for channel in PT_CHANNELS:
             plot(selected[selected["channel"] == channel], channel, selection, lumi,
                 os.path.join(args.outdir, f"zh_mass_{args.model}_{selection}_Zhpt{channel}.pdf"),
                 score_label(edges, args.years, channel, selection))


if __name__ == "__main__":
    main()
