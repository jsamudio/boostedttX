#!/usr/bin/env python3
"""
Summary figure of the DCTR hdamp gen-level check (matplotlib): does the TOP DCTR reweighting reproduce
the hdamp samples in the 5FS ttbar it was trained on, and in the 4FS ttbb?

Reads the .npz files written by dctr_gen_check.py (one per sample, twiki convention p0/p1 only).
One figure per observable (pT(tt), top pT): rows = samples, columns = hdamp up / down; hdamp sample /
nominal (points, stat.) and DCTR-reweighted nominal / nominal (line, stat. band), both unit area,
with chi2/ndf and the inclusive mean DCTR weight in each panel.

    python3 dctr_gen_check.py --sample TTToSemiLeptonic --era 2018
    python3 dctr_gen_check.py --sample TTbb_SemiLeptonic --era 2018
    python3 plot_dctr_summary.py plots/hdamp_dctr_gen/TTToSemiLeptonic_2018_dctr.npz \
        plots/hdamp_dctr_gen/TTbb_SemiLeptonic_2018_dctr.npz
    # tt+B: 5FS tt+B part vs 4FS (dctr_gen_check.py --ttb), e.g.
    python3 plot_dctr_summary.py plots/hdamp_dctr_gen/TTToSemiLeptonic_2018_ttB_dctr.npz \
        plots/hdamp_dctr_gen/TTbb_SemiLeptonic_2018_ttB_dctr.npz --tag ttB
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
       "pt_bb": r"$p_T(b\bar{b})$, additional [GeV]", "m_bb": r"$m(b\bar{b})$, additional [GeV]",
       "dr_bb": r"$\Delta R(b\bar{b})$, additional", "pt_ttbb": r"$p_T(t\bar{t}b\bar{b})$ [GeV]"}
LABEL = {"TTToSemiLeptonic": r"5FS $t\bar{t}$ (training process)", "TTbb_SemiLeptonic": r"4FS $t\bar{t}b\bar{b}$",
         "TTToSemiLeptonic tt+B": r"5FS $t\bar{t}$, tt+B part", "TTbb_SemiLeptonic tt+B": r"4FS $t\bar{t}b\bar{b}$, tt+B"}


def draw(runs, obs, outdir, ylim, tag=""):
    fig, axes = plt.subplots(len(runs), 2, figsize=(10, 3.4 * len(runs)), sharex=True, squeeze=False)
    for i, run in enumerate(runs):
        sample = str(run["sample"])
        for j, var in enumerate(("up", "down")):
            ax = axes[i, j]
            bins, rs, es, rd, ed = (run[f"{var}_{obs}_{k}"] for k in ("bins", "rs", "es", "rd", "ed"))
            x = 0.5 * (bins[1:] + bins[:-1])
            ax.errorbar(x, rs, yerr=es, fmt="ko", ms=4, label=f"hdamp {var} sample / nominal")
            ax.stairs(rd, bins, baseline=None, color="#d62728", lw=1.8, label="DCTR-reweighted nominal / nominal")
            ax.fill_between(bins, np.append(rd - ed, rd[-1] - ed[-1]), np.append(rd + ed, rd[-1] + ed[-1]),
                            step="post", color="#d62728", alpha=0.2, lw=0)
            ax.axhline(1, color="gray", lw=0.8)
            ax.set_ylim(*ylim)
            ax.text(0.03, 0.95, f"{LABEL.get(sample, sample)}, hdamp {var}\n"
                    rf"$\chi^2$/ndf = {float(run[f'{var}_{obs}_chi2']):.0f}/{int(run[f'{var}_{obs}_ndf'])}, "
                    f"mean DCTR weight {float(run[f'{var}_meanw']):.3f}"
                    + (f"\ntt+B fraction / nominal {float(run[f'{var}_ttbfrac_ratio']):.3f}"
                       if f"{var}_ttbfrac_ratio" in run else ""),
                    transform=ax.transAxes, va="top", fontsize=8)
            if j == 0:
                ax.set_ylabel("ratio to nominal (unit area)")
            if i == len(runs) - 1:
                ax.set_xlabel(OBS.get(obs, obs))
            if i == 0 and j == 0:
                ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"DCTR hdamp reweighting (TOP 13 TeV model), generator level, {str(runs[0]['era'])}", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{outdir}/dctr_summary_{obs}{'_' + tag if tag else ''}.{ext}")
    plt.close(fig)
    print(f"Wrote {outdir}/dctr_summary_{obs}{'_' + tag if tag else ''}.png/.pdf")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("npz", nargs="+", help="dctr_gen_check.py outputs, one per sample (rows in this order)")
    ap.add_argument("--obs", nargs="+", default=["pt_tt", "pt_top"], choices=list(OBS))
    ap.add_argument("--ylim", nargs=2, type=float, default=[0.8, 1.25])
    ap.add_argument("--outdir", default="plots/hdamp_dctr_gen")
    ap.add_argument("--tag", default="", help="suffix of the output names (dctr_summary_<obs>_<tag>)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    runs = [dict(np.load(f, allow_pickle=True)) for f in args.npz]
    for obs in args.obs:
        draw(runs, obs, args.outdir, args.ylim, args.tag)
