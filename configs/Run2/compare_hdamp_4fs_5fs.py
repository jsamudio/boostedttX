#!/usr/bin/env python3
"""
hdamp on tt+B: 4FS ttbb samples vs the tt+B part of the 5FS ttbar samples (Run 2, one era).

Question: does the tt+B part (genTtbarId % 100 in 51-56) of the 5FS ttbar hdamp samples reproduce the
hdamp effect of the dedicated 4FS ttbb hdamp samples in this selection? If so, the 2024 5FS ttbar hdamp
samples can stand in for the missing 2024 ttbb hdamp samples (the DCTR reweighting cannot:
review-notes "DCTR hdamp validation (ttbb 2018)").

Inputs (.coffea, analysis selection):
  --nominal-4fs  nominal TTbb_* (e.g. outputs_dctr/nominal/bkg)          subsamples TTbb_*__tt+B
  --hdamp-4fs    TTbb_*_hdampUP/DOWN (e.g. outputs_hdamp/.../variations) subsamples TTbb_*_hdamp*__tt+B
  --nominal-5fs  nominal TTTo* (e.g. outputs/nominal/bkg)                subsamples TTTo*__tt+B
  --hdamp-5fs    TTTo*_hdampUP/DOWN (variations group)                   subsamples TTTo*_hdamp*__tt+B
Each sample is normalized to the cross section of its nominal sample (outputToPkl.get_norm), so the
ratio variation / nominal is the hdamp effect on the selected tt+B yield (acceptance and, for 5FS, the
tt+B fraction). Selection makeDatacard.cuts, weight getZhbbWeight (compare_hdamp_dctr.read).

Output: rate ratio to nominal (all, Zhpt1-3) for 4FS and 5FS with stat. uncertainties and their
difference in sigma; per observable the shape ratio to nominal (unit area) of both and a chi2 of
4FS vs 5FS (independent samples). Plots <outdir>/<obs>.png/.pdf

    python3 compare_hdamp_4fs_5fs.py --era 2018 \
        --nominal-4fs outputs_dctr/nominal/bkg/output_all.coffea \
        --hdamp-4fs outputs_hdamp/nominal/variations/output_all.coffea \
        --nominal-5fs outputs/nominal/bkg/output_all.coffea \
        --hdamp-5fs outputs_hdamp_5fs/nominal/variations/output_all.coffea
"""
import argparse
import os

import numpy as np

from compare_hdamp_dctr import read, hist, ratio, yields, CHANNELS, OBSERVABLES, VARS

TTBB = ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu"]
TTBAR = ["TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu"]
SCHEMES = {"4FS": TTBB, "5FS": TTBAR}


def load(args):
    files = {"4FS": (args.nominal_4fs, args.hdamp_4fs), "5FS": (args.nominal_5fs, args.hdamp_5fs)}
    out = {}
    for scheme, bases in SCHEMES.items():
        nom_files, var_files = files[scheme]
        out[scheme] = {"nominal": read(nom_files, [f"{s}__tt+B" for s in bases], args.btag_norm, era=args.era)}
        for v, tag in VARS.items():
            out[scheme][v] = read(var_files, [f"{s}_{tag}__tt+B" for s in bases], args.btag_norm, era=args.era)
    return out


def channel_mask(df, lo, hi):
    return None if lo is None else ((df["ZH_pt"] > lo) & (df["ZH_pt"] <= hi)).to_numpy()


def rates(d):
    print(f"\n  {'region':7}" + "".join(f" {s + ' ' + VARS[v] + '/nom':>22}" for v in VARS for s in SCHEMES)
          + "".join(f" {'diff ' + v:>10}" for v in VARS))
    regions = {"all": (None, None), **CHANNELS}
    for name, (lo, hi) in regions.items():
        r, line = {}, f"  {name:7}"
        for v in VARS:
            for s in SCHEMES:
                nom, var = d[s]["nominal"], d[s][v]
                yn, en = yields(nom, nom["w"].to_numpy(), channel_mask(nom, lo, hi))
                yv, ev = yields(var, var["w"].to_numpy(), channel_mask(var, lo, hi))
                rr, re = ratio(np.array([yv]), np.array([ev]), np.array([yn]), np.array([en]))
                r[(s, v)] = (rr[0], re[0])
                line += f" {rr[0]:13.4f} +- {re[0]:.4f}"
        for v in VARS:
            (a, ea), (b, eb) = r[("4FS", v)], r[("5FS", v)]
            line += f" {(a - b) / np.hypot(ea, eb):9.1f}s"
        print(line)


def shapes(d, outdir, plot):
    print("\n  shape chi2 / ndf (4FS vs 5FS ratio to nominal, unit area):")
    for obs, (bins, label) in OBSERVABLES.items():
        if any(obs not in d[s][k] for s in SCHEMES for k in d[s]):
            continue
        curves, chi = {}, []
        for v in VARS:
            rr = {}
            for s in SCHEMES:
                hn, en = hist(d[s]["nominal"][obs], d[s]["nominal"]["w"].to_numpy(), bins)
                hv, ev = hist(d[s][v][obs], d[s][v]["w"].to_numpy(), bins)
                kn, kv = 1 / hn.sum(), 1 / hv.sum()
                rr[s] = ratio(hv * kv, ev * kv, hn * kn, en * kn)
            (a, ea), (b, eb) = rr["4FS"], rr["5FS"]
            ok = (a > 0) & (b > 0)
            c2 = (((a - b) ** 2) / np.maximum(ea ** 2 + eb ** 2, 1e-12))[ok].sum()
            chi.append(f"{v} {c2:.1f}/{max(int(ok.sum()) - 1, 0)}")
            curves[v] = rr
        print(f"    {obs:28} " + "  ".join(chi))
        if plot:
            draw(bins, label, curves, f"{outdir}/{obs}")


def draw(bins, label, curves, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    x = 0.5 * (bins[1:] + bins[:-1])
    dx = 0.12 * np.diff(bins)
    colors = {"up": "#d62728", "down": "#1f77b4"}
    vals = []
    for v, rr in curves.items():
        (a, ea), (b, eb) = rr["4FS"], rr["5FS"]
        ax.errorbar(x - dx, a, yerr=ea, fmt="o", ms=4, color=colors[v], label=f"4FS ttbb {VARS[v]} / nominal")
        ax.errorbar(x + dx, b, yerr=eb, fmt="s", ms=4, mfc="none", color=colors[v],
                    label=f"5FS tt+B {VARS[v]} / nominal")
        vals += list(a[a > 0]) + list(b[b > 0])
    span = min(max(0.1, np.abs(np.array(vals) - 1).max() if vals else 0.1), 0.6)
    ax.set_ylim(1 - 1.15 * span, 1 + 1.6 * span)
    ax.axhline(1, color="gray", lw=0.8)
    ax.set_xlabel(label)
    ax.set_ylabel("ratio to nominal (shape)")
    ax.legend(fontsize=7, ncol=2, loc="upper center")
    ax.set_title("tt+B hdamp: 4FS ttbb vs 5FS ttbar", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--era", default="2018")
    ap.add_argument("--nominal-4fs", nargs="+", required=True)
    ap.add_argument("--hdamp-4fs", nargs="+", required=True)
    ap.add_argument("--nominal-5fs", nargs="+", required=True)
    ap.add_argument("--hdamp-5fs", nargs="+", required=True)
    ap.add_argument("--btag-norm", default=None, help="see compare_hdamp_dctr.py")
    ap.add_argument("--outdir", default="plots/hdamp_4fs_5fs")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    d = load(args)
    for s in SCHEMES:
        print(f"{s} tt+B selected events: nominal {len(d[s]['nominal'])}, "
              + ", ".join(f"{VARS[v]} {len(d[s][v])}" for v in VARS))
    rates(d)
    shapes(d, args.outdir, not args.no_plots)
