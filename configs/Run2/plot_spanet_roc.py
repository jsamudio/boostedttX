#!/usr/bin/env python3
"""
ROC curves and normalized score shapes of a SPANet model on the analysis samples (Run 2).

Reads the analysis pickles (outputToPkl.py: PICKLE_DIR / <process>.pkl) and uses the discriminant
disc_<model> = sig / (sig + ttbb + ttlf), the event weight of the datacards (getZhbbWeight) and,
by default, the datacard selection (makeDatacard.cuts: >= 5 AK4, 2 b outside the candidate,
Xbb >= LP, MET, 50 < ZH_M < 200).

    python3 plot_spanet_roc.py                          # active model, all years + combined
    python3 plot_spanet_roc.py --model noOverlap_MultiWP --years 2018 --selection none

Outputs (in --outdir): roc_<model>_<year>.pdf/.png, score_<model>_<year>.pdf/.png, and the AUCs on
stdout. Signal = ttH + ttZ; ROC backgrounds: all, tt+B, tt+LF/C. CMS style and labels from mplhep.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import pandas as pd

from analysis_years import YEARS, PICKLE_DIR
from datacard_shapes import getZhbbWeight
from spanet_models import ACTIVE_MODELS

# process groups of the shape plot, in a fixed color order (signal solid, background dashed)
GROUPS = {
    "ttH":        {"processes": ["ttH"], "label": r"$t\bar{t}H$", "signal": True},
    "ttZ":        {"processes": ["ttZ"], "label": r"$t\bar{t}Z$", "signal": True},
    "tt_B":       {"processes": ["tt_B"], "label": r"$t\bar{t}+B$", "signal": False},
    "TTBar":      {"processes": ["TTBar"], "label": r"$t\bar{t}$+LF/C", "signal": False},
    "minor":      {"processes": ["SingleTop", "ttX", "VJets", "QCD"], "label": "Minor (t, ttX, V+jets, QCD)",
                   "signal": False},
}
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SIGNAL = ["ttH", "ttZ"]
ROC_BACKGROUNDS = {"all": ["tt_B", "TTBar", "SingleTop", "ttX", "VJets", "QCD"],
                   "tt_B": ["tt_B"], "TTBar": ["TTBar"]}
ROC_LABELS = {"all": "all backgrounds", "tt_B": r"$t\bar{t}+B$", "TTBar": r"$t\bar{t}$+LF/C"}


def load(years, model, selection):
    """One DataFrame per year with process, score and weight."""
    if selection == "datacard":
        from makeDatacard import cuts
    score = f"disc_{model}"
    out = {}
    processes = sorted({p for g in GROUPS.values() for p in g["processes"]})
    for year in years:
        frames = []
        for proc in processes:
            path = os.path.join(PICKLE_DIR.format(year=year), f"{proc}.pkl")
            if not os.path.exists(path):
                continue
            df = pd.read_pickle(path)
            if score not in df:
                raise SystemExit(f"{path} has no column {score}: rerun outputToPkl.py with the model active")
            if selection == "datacard":
                df = df[cuts(df, year)]
            frames.append(pd.DataFrame({"process": proc, "score": df[score].to_numpy(),
                                        "weight": getZhbbWeight(df, year).to_numpy()}))
        if frames:
            out[year] = pd.concat(frames, ignore_index=True)
        else:
            print(f"  [skip] no pickles for {year} in {PICKLE_DIR.format(year=year)}")
    return out


def roc(score_s, w_s, score_b, w_b):
    """(signal efficiency, background efficiency, AUC) scanning the cut on the score."""
    cuts = np.unique(np.concatenate([score_s, score_b]))[::-1]
    def eff(score, w):
        order = np.argsort(score)[::-1]
        cum = np.concatenate([[0.], np.cumsum(w[order])]) / w.sum()
        idx = np.searchsorted(-score[order], -cuts, side="right")      # events with score >= cut
        return cum[idx]
    eff_s, eff_b = np.concatenate([[0.], eff(score_s, w_s)]), np.concatenate([[0.], eff(score_b, w_b)])
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    auc = float(trapezoid(eff_s, eff_b))
    return eff_s, eff_b, auc


def plot_year(df, year, model, outdir, bins):
    tag = f"{model}_{year}"
    sig = df[df.process.isin(SIGNAL)]
    aucs = {}

    # ---- ROC
    fig, ax = plt.subplots()
    for (name, procs), color, ls in zip(ROC_BACKGROUNDS.items(), ["#2a78d6", "#1baf7a", "#eda100"], ["-", "--", ":"]):
        bkg = df[df.process.isin(procs)]
        if bkg.empty or sig.empty:
            continue
        eff_s, eff_b, auc = roc(sig.score.to_numpy(), sig.weight.to_numpy(),
                                bkg.score.to_numpy(), bkg.weight.to_numpy())
        aucs[name] = auc
        ax.plot(eff_s, eff_b, color=color, ls=ls, lw=3,
                label=f"{ROC_LABELS[name]} (AUC = {auc:.3f})")
    ax.plot([0, 1], [0, 1], color="#b0b0a8", lw=1, ls="-", zorder=0)
    ax.set_xlabel(r"Signal efficiency ($t\bar{t}H + t\bar{t}Z$)")
    ax.set_ylabel("Background efficiency")
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_ylim(1e-3, 1.05)
    ax.legend(loc="lower right", title=f"SPANet {model}")
    _label(ax, year)
    _save(fig, outdir, f"roc_{tag}")

    # ---- normalized score shapes
    fig, ax = plt.subplots()
    edges = np.linspace(0, 1, bins + 1)
    for (name, g), color in zip(GROUPS.items(), COLORS):
        sel = df[df.process.isin(g["processes"])]
        if sel.empty or sel.weight.sum() == 0:
            continue
        h, _ = np.histogram(sel.score, edges, weights=sel.weight)
        h2, _ = np.histogram(sel.score, edges, weights=sel.weight ** 2)
        norm = h.sum()
        hep.histplot(h / norm, edges, yerr=np.sqrt(h2) / norm, histtype="step", color=color, lw=3,
                     ls="-" if g["signal"] else "--", label=g["label"], ax=ax)
    ax.set_xlabel("SPANet score")
    ax.set_ylabel("a.u.")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.5 * ax.get_ylim()[1])        # headroom for the legend
    ax.legend(loc="upper center", ncol=2, fontsize="small", title=f"SPANet {model}")
    _label(ax, year)
    _save(fig, outdir, f"score_{tag}")
    return aucs


def _label(ax, year):
    """'CMS Simulation Preliminary' + era / sqrt(s) (normalized shapes: no luminosity)."""
    hep.cms.label("Preliminary", data=False, rlabel=f"{'Run 2' if year == 'Run2' else year} (13 TeV)", ax=ax)


def _save(fig, outdir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight", dpi=150)
    print(f"  wrote {os.path.join(outdir, name)}.pdf/.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=ACTIVE_MODELS[0] if ACTIVE_MODELS else "noOverlap_MultiWP")
    ap.add_argument("--years", nargs="+", default=list(YEARS))
    ap.add_argument("--selection", choices=["datacard", "none"], default="datacard",
                    help="datacard: makeDatacard.cuts (incl. Xbb >= LP); none: all pickled events (btag_mask)")
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--outdir", default="plots/spanet")
    ap.add_argument("--no-combined", action="store_true", help="skip the plots of all years together")
    args = ap.parse_args()

    plt.style.use(hep.style.CMS)

    os.makedirs(args.outdir, exist_ok=True)
    data = load(args.years, args.model, args.selection)
    if not data:
        raise SystemExit("No pickles found")
    if not args.no_combined and len(data) > 1:
        data["Run2"] = pd.concat(data.values(), ignore_index=True)

    print(f"\nSPANet {args.model}, selection: {args.selection}")
    table = {}
    for year, df in data.items():
        print(f"{year}: " + ", ".join(f"{p} {n}" for p, n in df.process.value_counts().items()))
        table[year] = plot_year(df, year, args.model, args.outdir, args.bins)
    names = {"all": "all bkg", "tt_B": "tt+B", "TTBar": "tt+LF/C"}
    print("\nAUC, signal = ttH + ttZ")
    print(f"{'':>8}" + "".join(f"{names[b]:>10}" for b in ROC_BACKGROUNDS))
    for year, aucs in table.items():
        print(f"{year:>8}" + "".join(f"{aucs.get(b, float('nan')):>10.4f}" for b in ROC_BACKGROUNDS))


if __name__ == "__main__":
    main()
