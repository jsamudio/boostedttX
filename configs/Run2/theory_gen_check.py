#!/usr/bin/env python3
"""
Generator level, 4FS ttbb: how pT(ttbb) and pT(bb) (additional b pair) change under the tt_B theory
uncertainties of the card, next to the hdamp samples. No selection (NanoAOD directly), tt+B events.

Theory variations (weights of the nominal sample, as in the analysis; dctr_gen_check.theory_weights):
  mu_r, mu_f   LHEScaleWeight                     isr, fsr   PSWeight
  pdf          LHEPdfWeight replicas, pdf_uncertainty.pdf_up_down per bin (as in the datacards)
  toppt        top pT reweighting: up = reweighting off, down = nominal (one-sided, as in the card)
The nominal weight is sign(genWeight) x top pT weight (the analysis nominal); hdamp samples likewise.

Normalization: each variation divided by its own sum of weights over ALL kept (tt+B) events, i.e.
the cross-section change is removed as in the card (process norms), and what remains is the change
of the fraction of events per bin (shape, and the fraction with >= 2 additional b jets where pT(bb)
and pT(ttbb) are defined). The hdamp samples are normalized the same way (same cross section).

Per observable: <obs>_theory.png with two panels (up, down): variation / nominal per bin, one line
per theory uncertainty (no stat. error: same events as the nominal) and the hdamp samples as points
with the stat. of both samples. Per-bin ratios printed.

    python3 theory_gen_check.py --sample TTbb_SemiLeptonic --era 2018 --max-events 1000000
"""
import argparse
import os

import numpy as np

from dctr_gen_check import read, files_of, observables, OBS
from pdf_uncertainty import pdf_up_down

THEORY = ["mu_r", "mu_f", "isr", "fsr", "pdf", "toppt"]


def hist(x, w, bins, total):
    """Histogram of the events where x is defined, divided by `total` (sum of weights of all events)."""
    fin = np.isfinite(x)
    xc = np.clip(x[fin], bins[0], bins[-1] - 1e-9)
    h, _ = np.histogram(xc, bins, weights=w[fin])
    e2, _ = np.histogram(xc, bins, weights=w[fin] ** 2)
    return h / total, np.sqrt(e2) / total


def theory_ratios(x, ev, bins):
    """{(name, 'up'/'down'): variation / nominal per bin}, each normalized to its own total."""
    w0 = ev["w"] * ev["th_toppt"]                        # analysis nominal
    hn, _ = hist(x, w0, bins, w0.sum())
    out = {}
    for name in ("mu_r", "mu_f", "isr", "fsr"):
        for ud in ("up", "down"):
            w = w0 * ev[f"th_{name}_{ud}"]
            out[(name, ud)] = hist(x, w, bins, w.sum())[0] / hn
    w_off = ev["w"]                                      # top pT: up = off, down = nominal
    out[("toppt", "up")] = hist(x, w_off, bins, w_off.sum())[0] / hn
    out[("toppt", "down")] = np.ones_like(hn)
    pdf = ev["th_pdf"]
    members = np.array([hist(x, w0 * pdf[:, k], bins, (w0 * pdf[:, k]).sum())[0] for k in range(pdf.shape[1])])
    up, down = pdf_up_down(members)
    out[("pdf", "up")], out[("pdf", "down")] = up / hn, down / hn
    return out, hn


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", default="TTbb_SemiLeptonic")
    ap.add_argument("--era", default="2018")
    ap.add_argument("--max-events", type=int, default=1000000)
    ap.add_argument("--obs", nargs="+", default=["pt_ttbb", "pt_bb"], choices=list(OBS))
    ap.add_argument("--outdir", default="plots/theory_gen")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    evs = {}
    for tag, name in (("nominal", args.sample), ("up", f"{args.sample}_hdampUP"), ("down", f"{args.sample}_hdampDOWN")):
        print(f"Reading {name} {args.era}")
        evs[tag], _ = read(files_of(name, args.era), args.max_events, ttb=True, theory=True)
    obs = {k: observables(v) for k, v in evs.items()}
    w_an = {k: v["w"] * v["th_toppt"] for k, v in evs.items()}      # analysis nominal weight of each sample

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = dict(zip(THEORY, plt.get_cmap("tab10").colors))

    for o in args.obs:
        bins, label = OBS[o]
        x0 = obs["nominal"][o]
        th, hn = theory_ratios(x0, evs["nominal"], bins)
        _, en = hist(x0, w_an["nominal"], bins, w_an["nominal"].sum())
        hd = {}
        for v in ("up", "down"):           # hdamp sample / nominal, same weight and normalization
            hv, ev_ = hist(obs[v][o], w_an[v], bins, w_an[v].sum())
            r = np.divide(hv, hn, out=np.full_like(hv, np.nan), where=hn > 0)
            er = r * np.sqrt(np.divide(ev_, hv, out=np.zeros_like(hv), where=hv > 0) ** 2
                             + np.divide(en, hn, out=np.zeros_like(hn), where=hn > 0) ** 2)
            hd[v] = (r, er)

        print(f"\n=== {o} (per bin, variation / nominal; fraction of all tt+B events with the observable "
              f"defined: {hn.sum():.4f})")
        print("    bins: " + " ".join(f"{b:g}" for b in bins))
        for v in ("up", "down"):
            for name in THEORY:
                print(f"    {name:6} {v:4}: " + " ".join(f"{r:.3f}" for r in th[(name, v)]))
            r, er = hd[v]
            print(f"    hdamp  {v:4}: " + " ".join(f"{a:.3f}({b:.3f})" for a, b in zip(r, er)))

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
        x = 0.5 * (bins[1:] + bins[:-1])
        for ax, v in zip(axes, ("up", "down")):
            for name in THEORY:
                ax.stairs(th[(name, v)], bins, baseline=None, color=colors[name], lw=1.5, label=name)
            r, er = hd[v]
            ax.errorbar(x, r, yerr=er, fmt="ko", ms=4, label="hdamp sample (stat.)")
            ax.axhline(1, color="gray", lw=0.8)
            ax.set_xlabel(label)
            ax.set_title(f"{v} variations", fontsize=10)
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("variation / nominal")
        axes[0].set_ylim(0.7, 1.3)
        axes[0].legend(fontsize=7, ncol=2, frameon=False)
        fig.suptitle(f"{args.sample} {args.era}, tt+B, generator level: theory uncertainties and hdamp "
                     f"(each normalized to its own total)", fontsize=9)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{args.outdir}/{args.sample}_{args.era}_{o}_theory.{ext}")
        plt.close(fig)
    print(f"\nWrote {args.outdir}/{args.sample}_{args.era}_<obs>_theory.png/.pdf")


if __name__ == "__main__":
    main()
