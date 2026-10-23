#!/usr/bin/env python3
"""
Pre-fit datacard templates (MC only, blinded): the flattened 2D (SPANet bin x Z/H mass bin)
templates of every channel side by side, stacked processes, signal on top, and the pre-fit
uncertainty band. One figure per datacard (year). Works for any datacard written by
makeDatacard.py (Run 2, and Run 3 once it has datacards in the same format).

    cd configs/Run2
    python3 ../common/plot_templates.py datacards/datacard_nominal_disc_<model>_2018.txt \\
        --edges nn_edges_disc_<model>.json --year 2018 --lumi 59.8 -o plots/templates

The uncertainty band per bin is the quadratic sum of
  * the MC statistical uncertainty (sum of the template variances),
  * every lnN and shape nuisance of the card, each fully correlated across processes:
    shift = sum over processes of (up - nominal) and (down - nominal), symmetrized,
excluding rateParams (e.g. the free tt+bb normalization). --stat-only shows only the first.
This approximates the combine pre-fit uncertainty, without correlations between bins.
"""
import argparse
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot

# card process -> plot group (regular expressions, first match), stacked bottom to top
GROUPS = [
    ("TTBar",     r"^TTBar$",          r"$t\bar{t}$+LF, $t\bar{t}+c\bar{c}$", "#ff9933"),
    ("tt_B",      r"^tt_B$",           r"$t\bar{t}+b\bar{b}$",                 "#4daf4a"),
    ("SingleTop", r"^SingleTop$",      "single top",                           "#8c564b"),
    ("ttX",       r"^ttX$",            r"$t\bar{t}t\bar{t}$, $t\bar{t}\gamma$, $t\bar{t}W$", "#b2182b"),
    ("VJets",     r"^VJets$",          "V+jets",                               "#66ccee"),
    ("QCD",       r"^QCD$",            "QCD",                                  "#999999"),
    ("ttH",       r"^ttH\d*$",         r"$t\bar{t}H$",                         "#e600e6"),
    ("ttZ",       r"^ttZ\d*$",         r"$t\bar{t}Z$",                         "#1f3a93"),
]
# Run 2 channels (makeDatacard.pt_bins = [0, 200, 300, 450])
CHANNEL_LABELS = {"Zhpt1": r"$200 < p_{T}^{Z/H\ cand.} < 300$",
                  "Zhpt2": r"$300 < p_{T}^{Z/H\ cand.} < 450$",
                  "Zhpt3": r"$p_{T}^{Z/H\ cand.} > 450$"}


# ---------------------------------------------------------------- datacard
def parse_card(path):
    """(columns, nuisances, shapes file): columns = [(channel, process)], nuisances =
    [(name, type, [value per column])] for lnN and shape nuisances."""
    lines = [l.split("#")[0].strip() for l in open(path)]
    lines = [l for l in lines if l and not set(l) <= {"-"}]
    shapes_file, bins = None, None
    columns, nuisances = None, []
    for i, line in enumerate(lines):
        tok = line.split()
        if tok[0] == "shapes" and tok[1] == "*" and shapes_file is None:
            shapes_file = os.path.join(os.path.dirname(path), tok[3])
        elif tok[0] == "bin" and i + 1 < len(lines) and lines[i + 1].split()[0] == "process":
            bins = tok[1:]
            columns = list(zip(bins, lines[i + 1].split()[1:]))
        elif columns and len(tok) >= 2 + len(columns) and tok[1] in ("lnN", "shape"):
            nuisances.append((tok[0], tok[1], tok[2:2 + len(columns)]))
    if columns is None:
        raise SystemExit(f"No bin/process lines found in {path}")
    return columns, nuisances, shapes_file


def lnn(value):
    """(kappa_down, kappa_up) of an lnN entry ('1.05' or '0.95/1.06'); None for '-'."""
    if value == "-":
        return None
    if "/" in value:
        down, up = (float(v) for v in value.split("/"))
        return down, up
    k = float(value)
    return 1 / k, k


# ---------------------------------------------------------------- templates
def load(card):
    """{channel: {"groups": {group: sumw}, "stat2": array, "shifts": [array per nuisance]}}."""
    columns, nuisances, shapes_file = parse_card(card)
    f = uproot.open(shapes_file)
    channels = list(dict.fromkeys(ch for ch, _ in columns))
    out = {}
    for ch in channels:
        nominal, groups = {}, {}
        stat2 = None
        for (c, proc) in columns:
            if c != ch:
                continue
            h = f[f"{ch}_{proc}"]
            nominal[proc] = h.values()
            var = h.variances()
            stat2 = var if stat2 is None else stat2 + var
            group = next((g for g, pattern, _, _ in GROUPS if re.match(pattern, proc)), "other")
            groups[group] = groups.get(group, 0) + nominal[proc]
        shifts = []
        for name, kind, values in nuisances:
            up = np.zeros_like(stat2)
            down = np.zeros_like(stat2)
            for (c, proc), value in zip(columns, values):
                if c != ch or value == "-":
                    continue
                if kind == "lnN":
                    kd, ku = lnn(value)
                    up += nominal[proc] * (ku - 1)
                    down += nominal[proc] * (kd - 1)
                else:
                    key = f"{ch}_{proc}_{name}"
                    if f"{key}Up" in f and f"{key}Down" in f:
                        up += f[f"{key}Up"].values() - nominal[proc]
                        down += f[f"{key}Down"].values() - nominal[proc]
            shifts.append(0.5 * (np.abs(up) + np.abs(down)))
        out[ch] = {"groups": groups, "stat2": stat2, "shifts": shifts}
    return out


def nn_bins(edges_file, year, channel, n_bins):
    """Number of SPANet bins of a channel (from the edges JSON; 6 if not given)."""
    if edges_file:
        edges = json.load(open(edges_file))
        e = edges.get(year) or next(iter(edges.values()))
        return len(e[channel]) - 1
    return 6 if n_bins % 6 == 0 else n_bins


# ---------------------------------------------------------------- plot
def plot(templates, out, edges_file, year, lumi, stat_only, com):
    plt.style.use(hep.style.CMS)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(18, 9), gridspec_kw={"height_ratios": [4, 1]}, sharex=True)
    channels = list(templates)
    present = [g for g, *_ in GROUPS if any(g in templates[ch]["groups"] for ch in channels)]
    others = any("other" in templates[ch]["groups"] for ch in channels)

    stacks = {g: [] for g in present + (["other"] if others else [])}
    total, unc = [], []
    boundaries, nn_lines, labels = [0], [], []
    for ch in channels:
        t = templates[ch]
        n = len(t["stat2"])
        for g in stacks:
            stacks[g].append(t["groups"].get(g, np.zeros(n)))
        tot = sum(t["groups"].values())
        var = t["stat2"] + (0 if stat_only else sum(s ** 2 for s in t["shifts"]))
        total.append(tot)
        unc.append(np.sqrt(var))
        n_nn = nn_bins(edges_file, year, ch, n)
        n_mass = n // n_nn
        start = boundaries[-1]
        nn_lines += [start + k * n_mass for k in range(1, n_nn)]
        labels.append((start, start + n, ch, n_nn, n_mass))
        boundaries.append(start + n)
    total, unc = np.concatenate(total), np.concatenate(unc)
    x = np.arange(len(total) + 1)

    group_style = {g: (label, color) for g, _, label, color in GROUPS}
    group_style["other"] = ("other", "#cccccc")
    order = list(stacks)
    hep.histplot([np.concatenate(stacks[g]) for g in order], x, stack=True, histtype="fill",
                 label=[group_style[g][0] for g in order], color=[group_style[g][1] for g in order], ax=ax)
    band = "stat." if stat_only else "stat.+syst."
    ax.stairs(total + unc, x, baseline=np.clip(total - unc, 1e-6, None), fill=True, hatch="xxxx",
              facecolor="none", edgecolor="#d62728", linewidth=0, label=band)

    rel = np.divide(unc, total, out=np.zeros_like(unc), where=total > 0)
    rax.stairs(1 + rel, x, baseline=1 - rel, fill=True, hatch="xxxx", facecolor="none",
               edgecolor="#d62728", linewidth=0)
    rax.axhline(1, color="black", lw=1)
    rax.set_ylim(0.5, 1.5)
    rax.set_ylabel("Unc. / MC", fontsize="small")
    rax.set_xlabel("Analysis bins")

    ax.set_yscale("log")
    ymax = np.max(total) if np.any(total > 0) else 1
    ax.set_ylim(max(1e-2, 1e-3 * ymax), 1e3 * ymax)
    ax.set_ylabel("Events / bin")
    for a in (ax, rax):
        for b in boundaries[1:-1]:
            a.axvline(b, color="black", lw=1.5)
        for b in nn_lines:
            a.axvline(b, color="black", lw=0.8, ls="--")
    for start, stop, ch, n_nn, n_mass in labels:
        ax.text(0.5 * (start + stop), 3e2 * ymax, CHANNEL_LABELS.get(ch, ch), ha="center", fontsize="x-small")
        if start == 0:
            ax.text(start + 0.3, 30 * ymax, r"NN$_{bin}$", fontsize="x-small")
            for k in range(n_nn):
                ax.text(start + (k + 0.5) * n_mass, 8 * ymax, str(k + 1), ha="center", fontsize="x-small")
    ax.set_xlim(0, len(total))
    ax.legend(ncol=3, fontsize="x-small", loc="upper right", bbox_to_anchor=(1.0, 0.86))
    hep.cms.label("Preliminary", data=False, lumi=lumi, com=com, ax=ax)

    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cards", nargs="+", help="datacard .txt file(s), one figure each")
    ap.add_argument("--edges", default=None, help="nn_edges_disc_<model>.json (number of SPANet bins)")
    ap.add_argument("--year", default=None, help="year key in the edges file (default: from the card name)")
    ap.add_argument("--lumi", type=float, default=None, help="luminosity label [fb^-1]")
    ap.add_argument("--com", type=float, default=13, help="centre-of-mass energy label [TeV]")
    ap.add_argument("--stat-only", action="store_true", help="band = MC statistical uncertainty only")
    ap.add_argument("-o", "--outdir", default="plots/templates")
    args = ap.parse_args()
    for card in args.cards:
        year = args.year or (re.findall(r"(20\d\d)", os.path.basename(card)) or [None])[-1]
        out = os.path.join(args.outdir, os.path.splitext(os.path.basename(card))[0] + ".pdf")
        plot(load(card), out, args.edges, year, args.lumi, args.stat_only, args.com)


if __name__ == "__main__":
    main()
