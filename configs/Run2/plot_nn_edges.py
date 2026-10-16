#!/usr/bin/env python3
"""
Plot the SPANet (NN) bin edges of the datacards: one row per year and ZH-candidate pT channel,
one marker per edge. The edges come from DataCardShapes.save_edges (makeDatacard.py writes
nn_edges_<disc>.json and calls plot_edges on it), or any JSON {year: {"Zhpt1": [...], ...}}.

    python3 plot_nn_edges.py nn_edges_disc_noOverlap_MultiWP.json [-o plots/nn_edges.pdf]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

# reco ZH-candidate pT channels of the datacards (makeDatacard.pt_bins = [0, 200, 300, 450])
PT_CHANNELS = {"Zhpt1": r"$200 < p_{T}^{ZH\ cand.} < 300$",
               "Zhpt2": r"$300 < p_{T}^{ZH\ cand.} < 450$",
               "Zhpt3": r"$p_{T}^{ZH\ cand.} > 450$"}
COLORS = {"Zhpt1": "#2a78d6", "Zhpt2": "#eb6834", "Zhpt3": "#1baf7a"}


def plot_edges(edges, out, title="Analysis quantile-defined SPANet bins", xlabel="SPANet bin edge"):
    """edges: {year: {channel: [edge, ...]}}; writes `out` (.pdf) and the .png next to it."""
    plt.style.use(hep.style.CMS)
    rows = [(year, ch) for year in sorted(edges) for ch in sorted(edges[year])]
    fig, ax = plt.subplots(figsize=(12, 1.0 + 0.75 * len(rows)))
    for i, (year, ch) in enumerate(rows):
        e = np.asarray(edges[year][ch], dtype=float)
        ax.plot(e, np.full(len(e), i), ls="none", marker="*", ms=16, color=COLORS.get(ch, "#4a3aa7"),
                markeredgecolor="white", markeredgewidth=0.8, zorder=3)
        ax.axhline(i, color="#b0b0a8", lw=0.8, zorder=1)
    # separate the years
    years = [y for y, _ in rows]
    for i in range(1, len(rows)):
        if years[i] != years[i - 1]:
            ax.axhline(i - 0.5, color="black", lw=1.2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{y}  {PT_CHANNELS.get(ch, ch)}" for y, ch in rows], fontsize="x-small")
    ax.tick_params(axis="y", which="minor", left=False, right=False)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-0.02, 1.02)
    ax.xaxis.grid(True, color="#b0b0a8", lw=0.8)
    ax.set_xlabel(xlabel)
    hep.cms.label("Preliminary", data=False, ax=ax, rlabel="(13 TeV)")
    fig.suptitle(title, fontsize="medium", y=1.0 + 0.25 / fig.get_size_inches()[1])
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved NN bin edge plot to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("edges", help="JSON {year: {channel: [edges]}} (nn_edges_<disc>.json)")
    ap.add_argument("-o", "--out", default=None, help="output .pdf (default: plots/<json name>.pdf)")
    ap.add_argument("--title", default="Analysis quantile-defined SPANet bins")
    args = ap.parse_args()
    with open(args.edges) as f:
        edges = json.load(f)
    out = args.out or os.path.join("plots", os.path.splitext(os.path.basename(args.edges))[0] + ".pdf")
    plot_edges(edges, out, args.title)


if __name__ == "__main__":
    main()
