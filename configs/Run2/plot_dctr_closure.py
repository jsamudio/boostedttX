#!/usr/bin/env python3
"""
Closure figure of the hdamp reweighting at generator level (matplotlib): for each observable, the
distributions of the nominal sample, the dedicated hdamp up / down samples and the nominal reweighted
to up / down, all on one panel, with the ratio to the nominal below.

Reads one .npz written by dctr_gen_check.py (p0/p1 convention; 5FS model, or the 4FS ttbb model with
--model-4f, file name *_4f_dctr.npz). Up to four observables per figure, one per quadrant; each quadrant
= distribution (unit area / bin width; the last bin includes the overflow) + ratio to nominal.
Points: dedicated samples (stat.); dashed lines: reweighted nominal (stat. band in the ratio); grey band:
nominal stat. chi2/ndf of each
sample against its reweighted nominal (stat. of both) and the inclusive mean weight in the legend.

    python3 dctr_gen_check.py --sample TTbb_SemiLeptonic --era 2018 --model-4f <dir>
    python3 plot_dctr_closure.py plots/hdamp_dctr_gen/TTbb_SemiLeptonic_2018_4f_dctr.npz
    python3 plot_dctr_closure.py <npz> --obs pt_ttbb pt_ttbb_lhe pt_tt pt_bb --logy
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OBS = {"pt_tt": r"$p_T(t\bar{t})$ [GeV]", "pt_top": r"top $p_T$ [GeV]",
       "m_tt": r"$m(t\bar{t})$ [GeV]", "y_tt": r"$y(t\bar{t})$",
       "n_addb": "N additional b jets", "pt_b1": r"leading additional b jet $p_T$ [GeV]",
       "pt_bb": r"$p_T(b\bar{b})$, additional b jets [GeV]", "m_bb": r"$m(b\bar{b})$, additional b jets [GeV]",
       "dr_bb": r"$\Delta R(b\bar{b})$, additional b jets",
       "pt_ttbb": r"$p_T(t\bar{t}b\bar{b})$, additional b jets [GeV]",
       "pt_ttbb_lhe": r"$p_T(t\bar{t}b\bar{b})$, LHE b quarks [GeV]",
       "n_addj": "N additional light jets", "pt_j1": r"leading additional light jet $p_T$ [GeV]"}
LABEL = {"TTToSemiLeptonic": r"5FS $t\bar{t}$", "TTbb_SemiLeptonic": r"4FS $t\bar{t}b\bar{b}$"}
COLOR = {"up": "#d62728", "down": "#1f77b4"}
MARKER = {"up": "^", "down": "v"}


def quadrant(fig, spec, run, obs, logy, ylim):
    """Distribution (top) and ratio to nominal (bottom) of one observable in the given subplot spec."""
    gs = spec.subgridspec(2, 1, height_ratios=[3, 1.3], hspace=0.05)
    ax, rax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    rax.sharex(ax)
    bins = run[f"up_{obs}_bins"]
    width, x = np.diff(bins), 0.5 * (bins[1:] + bins[:-1])
    hn, en = run[f"nominal_{obs}_h"], run[f"nominal_{obs}_e"]
    ax.stairs(hn / width, bins, color="black", lw=1.6, baseline=None, label="nominal")
    rax.fill_between(bins, np.append(1 - en / hn, 1 - en[-1] / hn[-1]), np.append(1 + en / hn, 1 + en[-1] / hn[-1]),
                     step="post", color="gray", alpha=0.25, lw=0)
    for var in ("up", "down"):
        hs, es = run[f"{var}_{obs}_h"], run[f"{var}_{obs}_e"]
        rd, ed = run[f"{var}_{obs}_rd"], run[f"{var}_{obs}_ed"]
        chi2, ndf = float(run[f"{var}_{obs}_chi2"]), int(run[f"{var}_{obs}_ndf"])
        c = COLOR[var]
        ax.errorbar(x, hs / width, yerr=es / width, fmt=MARKER[var], color=c, ms=4, lw=1,
                    label=f"hdamp {var} sample")
        ax.stairs(rd * hn / width, bins, color=c, lw=1.6, ls="--", baseline=None,
                  label=rf"reweighted {var} ($\chi^2$/ndf {chi2:.0f}/{ndf}, $\langle w\rangle$ "
                        f"{float(run[f'{var}_meanw']):.3f})")
        rax.errorbar(x, run[f"{var}_{obs}_rs"], yerr=run[f"{var}_{obs}_es"], fmt=MARKER[var], color=c, ms=4, lw=1)
        rax.stairs(rd, bins, color=c, lw=1.6, ls="--", baseline=None)
        rax.fill_between(bins, np.append(rd - ed, rd[-1] - ed[-1]), np.append(rd + ed, rd[-1] + ed[-1]),
                         step="post", color=c, alpha=0.15, lw=0)
    rax.axhline(1, color="gray", lw=0.8)
    if logy:
        ax.set_yscale("log")
    ax.set_ylabel("unit area / bin width")
    ax.legend(fontsize=6.5, loc="upper right" if logy else "best")
    ax.tick_params(labelbottom=False)
    if ylim is None:                            # cover every ratio point and line, with a margin
        vals = np.concatenate([np.concatenate([run[f"{v}_{obs}_rs"] - run[f"{v}_{obs}_es"],
                                               run[f"{v}_{obs}_rs"] + run[f"{v}_{obs}_es"], run[f"{v}_{obs}_rd"]])
                               for v in ("up", "down")])
        vals = vals[np.isfinite(vals) & (vals > 0)]
        lo, hi = min(vals.min(), 0.95), max(vals.max(), 1.05)
        ylim = (max(lo - 0.08 * (hi - lo), 0), hi + 0.08 * (hi - lo))
    rax.set_ylim(*ylim)
    rax.set_ylabel("ratio to\nnominal", fontsize=8)
    rax.set_xlabel(OBS.get(obs, obs))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz", help="dctr_gen_check.py output (*_dctr.npz)")
    ap.add_argument("--obs", nargs="+", default=["pt_ttbb", "pt_ttbb_lhe", "pt_tt", "pt_bb"], choices=list(OBS),
                    help="observables, four per figure (quadrants); missing ones are skipped")
    ap.add_argument("--logy", action="store_true", help="log scale for the distributions")
    ap.add_argument("--ylim", nargs=2, type=float, help="ratio panel range (default: from the data)")
    ap.add_argument("--outdir", default="plots/hdamp_dctr_gen")
    ap.add_argument("--model-label", help="title label of the reweighting (default from the file name: 4FS "
                                          "ttbb classifier for *_4f_*, else TOP 5FS DCTR)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    run = dict(np.load(args.npz, allow_pickle=True))
    obs = [o for o in args.obs if f"up_{o}_bins" in run]
    for o in set(args.obs) - set(obs):
        print(f"skipped {o}: not in {args.npz}" + (" (needs --model-4f)" if o == "pt_ttbb_lhe" else ""))
    is4f = "_4f" in os.path.basename(args.npz)
    model = args.model_label or ("4FS ttbb classifier" if is4f else "TOP 5FS DCTR")
    sample = str(run["sample"])
    stem = os.path.basename(args.npz).replace("_dctr.npz", "")
    for k in range(0, len(obs), 4):
        group = obs[k:k + 4]
        ncol = 1 if len(group) == 1 else 2
        nrow = (len(group) + 1) // 2 if ncol == 2 else 1
        fig = plt.figure(figsize=(6.2 * ncol, 5.6 * nrow))
        outer = fig.add_gridspec(nrow, ncol, hspace=0.3, wspace=0.28)
        for i, o in enumerate(group):
            quadrant(fig, outer[i // ncol, i % ncol], run, o, args.logy, args.ylim)
        fig.suptitle(f"hdamp reweighting closure, generator level: {LABEL.get(sample, sample)} "
                     f"{str(run['era'])}, {model} (p0/p1)", fontsize=11)
        name = f"{args.outdir}/{stem}_closure_{'_'.join(group)}"
        for ext in ("png", "pdf"):
            fig.savefig(f"{name}.{ext}", bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {name}.png/.pdf")
