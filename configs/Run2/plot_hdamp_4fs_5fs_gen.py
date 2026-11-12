#!/usr/bin/env python3
"""
hdamp in the 4FS ttbb vs the 5FS ttbar (tt+B part), generator level, e.g. pT(ttbb).

Reads two .npz files written by dctr_gen_check.py --ttb (4FS and 5FS, same era), which hold unit-area
histograms of nominal, hdampUP and hdampDOWN per observable. Per observable:

  <obs>_nominal.png  nominal shape (unit area) of the 4FS and the 5FS tt+B, ratio 5FS / 4FS
  <obs>_hdamp.png    hdamp up and down: variation / nominal (unit area) of the 4FS and the 5FS tt+B

and per hdamp direction, printed and in the panels (statistical uncertainties of the variation and
the nominal sample, independent samples, bins with entries in all four histograms):
  flat   chi2 of variation / nominal against 1: is there any shape effect?
  slope  weighted linear fit of variation / nominal vs the bin centre, per 100 GeV, with its error:
         a trend test, more sensitive than chi2 to a monotonic hardening / softening
  4 vs 5 chi2 of the 4FS ratio against the 5FS ratio: is the hdamp shape effect the same in both?
  maxdev largest |ratio - 1| and its significance
Shapes only (unit area); the tt+B rate change is in the panel text (tt+B fraction / nominal).

    python3 dctr_gen_check.py --sample TTbb_SemiLeptonic --era 2018 --ttb --max-events 1000000
    python3 dctr_gen_check.py --sample TTToSemiLeptonic --era 2018 --ttb --max-events 5000000
    python3 plot_hdamp_4fs_5fs_gen.py plots/hdamp_dctr_gen/TTbb_SemiLeptonic_2018_ttB_dctr.npz \
        plots/hdamp_dctr_gen/TTToSemiLeptonic_2018_ttB_dctr.npz --obs pt_ttbb pt_bb pt_tt pt_j1
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from plot_dctr_summary import OBS  # noqa: E402

STYLE = {"4FS": dict(color="#1f4e9c", marker="o", mfc="#1f4e9c", label=r"4FS $t\bar{t}b\bar{b}$, tt+B"),
         "5FS": dict(color="#d9730d", marker="s", mfc="white", label=r"5FS $t\bar{t}$, tt+B part")}


def ratio(run, var, obs):
    """variation / nominal (unit area) and its error from both samples; mask of usable bins."""
    hn, en, hv, ev = (np.asarray(run[f"{t}_{obs}_{k}"], float)
                      for t, k in (("nominal", "h"), ("nominal", "e"), (var, "h"), (var, "e")))
    ok = (hn > 0) & (hv > 0)
    r = np.divide(hv, hn, out=np.ones_like(hn), where=ok)
    er = r * np.sqrt(np.divide(ev, hv, out=np.zeros_like(hv), where=ok) ** 2
                     + np.divide(en, hn, out=np.zeros_like(hn), where=ok) ** 2)
    return r, er, ok


def slope(x, r, er):
    """Weighted straight-line fit r = a + b (x - xm); b per 100 GeV and its error."""
    w = 1.0 / er ** 2
    xm = np.sum(w * x) / np.sum(w)
    sxx = np.sum(w * (x - xm) ** 2)
    b = np.sum(w * (x - xm) * r) / sxx
    return 100 * b, 100 / np.sqrt(sxx)


def tests(runs, var, obs, x):
    out, rs = {}, {}
    for name, run in runs.items():
        r, er, ok = ratio(run, var, obs)
        rs[name] = (r, er, ok)
        chi2 = np.sum(((r - 1) / er)[ok] ** 2)
        b, eb = slope(x[ok], r[ok], er[ok])
        dev = np.abs(r - 1)
        i = np.argmax(np.where(ok, dev, -1))
        out[name] = {"chi2": chi2, "ndf": int(ok.sum()) - 1, "slope": b, "eslope": eb,
                     "maxdev": dev[i], "maxsig": dev[i] / er[i], "maxbin": i}
    (r4, e4, k4), (r5, e5, k5) = rs["4FS"], rs["5FS"]
    both = k4 & k5
    out["4vs5"] = {"chi2": np.sum(((r4 - r5) ** 2 / (e4 ** 2 + e5 ** 2))[both]), "ndf": int(both.sum()) - 1}
    return out, rs


def draw_nominal(runs, obs, bins, outdir, era):
    x, wid = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(6, 5.5), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
    h = {}
    for name, run in runs.items():
        hn, en = (np.asarray(run[f"nominal_{obs}_{k}"], float) for k in ("h", "e"))
        h[name] = (hn, en)
        st = STYLE[name]
        ax.stairs(hn / wid, bins, color=st["color"], lw=1.6,
                  label=f"{st['label']} ({int(run[f'nominal_{obs}_n'])} events)")
        ax.errorbar(x, hn / wid, yerr=en / wid, fmt="none", ecolor=st["color"], lw=1)
    (h4, e4), (h5, e5) = h["4FS"], h["5FS"]
    ok = (h4 > 0) & (h5 > 0)
    r = np.divide(h5, h4, out=np.ones_like(h4), where=ok)
    er = r * np.sqrt(np.divide(e5, h5, out=np.zeros_like(h5), where=ok) ** 2
                     + np.divide(e4, h4, out=np.zeros_like(h4), where=ok) ** 2)
    chi2 = np.sum(((r - 1) / er)[ok] ** 2)
    axr.errorbar(x[ok], r[ok], yerr=er[ok], xerr=wid[ok] / 2, fmt="ko", ms=3.5, lw=1)
    axr.axhline(1, color="gray", lw=0.8)
    axr.set_ylabel("5FS / 4FS")
    axr.set_xlabel(OBS.get(obs, obs))
    ax.set_ylabel("fraction of events / GeV")
    ax.set_yscale("log")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(rf"Nominal shapes, generator level, {era}: $\chi^2$/ndf (5FS vs 4FS) = {chi2:.0f}/{int(ok.sum()) - 1}",
                 fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{outdir}/{obs}_nominal.{ext}")
    plt.close(fig)
    print(f"  nominal shape 5FS/4FS per bin: " + " ".join(f"{v:.3f}" for v in r)
          + f"   chi2 {chi2:.1f}/{int(ok.sum()) - 1}")


def draw_hdamp(runs, obs, bins, outdir, era, ylim):
    x, wid = 0.5 * (bins[1:] + bins[:-1]), np.diff(bins)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, var in zip(axes, ("up", "down")):
        res, rs = tests(runs, var, obs, x)
        print(f"  hdamp {var}")
        lines = []
        for k, name in enumerate(("4FS", "5FS")):
            r, er, ok = rs[name]
            st = STYLE[name]
            shift = (k - 0.5) * 0.12 * wid           # offset the points so the error bars do not overlap
            ax.errorbar(x[ok] + shift[ok], r[ok], yerr=er[ok], fmt=st["marker"], color=st["color"],
                        mfc=st["mfc"], ms=4.5, lw=1, label=st["label"])
            t = res[name]
            frac = runs[name].get(f"{var}_ttbfrac_ratio")
            lines.append(rf"{name}: flat $\chi^2$ {t['chi2']:.1f}/{t['ndf']}, slope {t['slope']:+.3f} $\pm$ "
                         f"{t['eslope']:.3f} /100 GeV" + (f", tt+B frac. {float(frac):.3f}" if frac is not None else ""))
            print(f"    {name} ratio per bin: " + " ".join(f"{v:.3f}" for v in r))
            print(f"    {name} error per bin: " + " ".join(f"{v:.3f}" for v in er))
            print(f"    {name} flat chi2 {t['chi2']:.1f}/{t['ndf']}, slope {t['slope']:+.4f} +- {t['eslope']:.4f} "
                  f"per 100 GeV ({t['slope'] / t['eslope']:+.1f} sigma), max |r-1| {t['maxdev']:.3f} "
                  f"({t['maxsig']:.1f} sigma, bin {bins[t['maxbin']]:g}-{bins[t['maxbin'] + 1]:g})")
        c = res["4vs5"]
        lines.append(rf"4FS vs 5FS: $\chi^2$ {c['chi2']:.1f}/{c['ndf']}")
        print(f"    4FS vs 5FS ratio chi2 {c['chi2']:.1f}/{c['ndf']}")
        ax.axhline(1, color="gray", lw=0.8)
        ax.set_ylim(*ylim)
        ax.set_xlabel(OBS.get(obs, obs))
        ax.set_title(f"hdamp {var}", fontsize=10)
        ax.text(0.03, 0.04, "\n".join(lines), transform=ax.transAxes, fontsize=7.5, va="bottom")
        if var == "up":
            ax.set_ylabel("hdamp variation / nominal (unit area)")
            ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle(f"hdamp shape effect, 4FS ttbb vs 5FS ttbar (tt+B), generator level, {era}", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{outdir}/{obs}_hdamp.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fs4", help="dctr_gen_check.py --ttb output of the 4FS ttbb")
    ap.add_argument("fs5", help="dctr_gen_check.py --ttb output of the 5FS ttbar")
    ap.add_argument("--obs", nargs="+", default=["pt_ttbb"], choices=list(OBS))
    ap.add_argument("--ylim", nargs=2, type=float, default=[0.75, 1.3])
    ap.add_argument("--outdir", default="plots/hdamp_4fs_5fs_gen")
    args = ap.parse_args()
    runs = {"4FS": dict(np.load(args.fs4, allow_pickle=True)), "5FS": dict(np.load(args.fs5, allow_pickle=True))}
    for name, run in runs.items():
        if "nominal_pt_tt_h" not in run:
            raise SystemExit(f"{name}: no nominal histograms in the npz; rerun dctr_gen_check.py (version with "
                             f"nominal_<obs>_h) on it")
        if "tt+B" not in str(run["sample"]):
            print(f"WARNING: {name} npz was not made with --ttb (all events, not only tt+B)")
    era = str(runs["4FS"]["era"])
    if str(runs["5FS"]["era"]) != era:
        raise SystemExit(f"different eras: 4FS {era}, 5FS {runs['5FS']['era']}")
    os.makedirs(args.outdir, exist_ok=True)
    for obs in args.obs:
        bins = np.asarray(runs["4FS"][f"up_{obs}_bins"], float)
        if not np.array_equal(bins, np.asarray(runs["5FS"][f"up_{obs}_bins"], float)):
            raise SystemExit(f"{obs}: different binning in the two files")
        print(f"=== {obs}")
        draw_nominal(runs, obs, bins, args.outdir, era)
        draw_hdamp(runs, obs, bins, args.outdir, era, args.ylim)
    print(f"Wrote {args.outdir}/<obs>_nominal.png, <obs>_hdamp.png (+ .pdf)")
