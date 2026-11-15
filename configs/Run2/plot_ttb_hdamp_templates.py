#!/usr/bin/env python3
"""
tt_B datacard templates: hdamp up / down (dedicated 4FS samples) next to the other tt_B theory
uncertainties of the card, bin by bin. One figure, one panel per pT channel:

  x       datacard bin (NN bin x mass bin, as in the card)
  y       variation / nominal, each variation scaled to the nominal total (shape only: the tt_B
          normalization floats in the fit, CMS_ttbbnorm), or with the rate kept (--keep-rate)
  band    the other tt_B theory uncertainties of the card in quadrature: per bin, for each shape
          nuisance the larger of |Up/Nom - 1| and |Down/Nom - 1|, summed in quadrature
  points  hdamp up and down / nominal from the samples, with the statistical uncertainty of the
          hdamp sample and the nominal (independent samples, so they scatter)

Reading: hdamp points inside the band (within their error bars) = hdamp shape change not larger than
what the card's other theory uncertainties already cover.

Inputs: the Run 2 card (.txt + .root, e.g. from sum_datacards.py) and the hdamp yields per datacard
bin from compare_hdamp_dctr.py --export (same model, years and NN edges). The hdamp Zhpt1 bins are
merged like the card (last two mass bins summed).

    python3 plot_ttb_hdamp_templates.py datacards/datacard_nominal_disc_<model>_Run2.txt \\
        plots/hdamp_run2/hdamp_ttB_<model>_<years>.json
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
from sum_datacards import parse_card, Templates  # noqa: E402
from compare_ttb_theory import card_nuisances  # noqa: E402


def to_card_bins(h, e2, n_nn, n_card):
    """hdamp bins (NN bin x 4 mass bins) -> card bins: identical, or the last two mass bins merged
    (Zhpt1 in the blinded card, makeDatacard.merge_last_mbin)."""
    h, e2 = np.asarray(h, float), np.asarray(e2, float)
    if len(h) == n_card:
        return h, e2
    n_mass = len(h) // n_nn
    if n_nn * (n_mass - 1) != n_card:
        raise SystemExit(f"cannot map {len(h)} hdamp bins onto {n_card} card bins")
    def merge(a):
        a = a.reshape(n_nn, n_mass).copy()
        a[:, -2] += a[:, -1]
        return a[:, :-1].ravel()
    return merge(h), merge(e2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("card")
    ap.add_argument("hdamp", help="compare_hdamp_dctr.py --export JSON")
    ap.add_argument("--process", default="tt_B")
    ap.add_argument("--keep-rate", action="store_true", help="do not scale the variations to the nominal total")
    ap.add_argument("--outdir", default="plots/ttb_theory")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    card = parse_card(args.card)
    tpl = Templates(card)
    channels = list(dict.fromkeys(ch for ch, p, _ in card["columns"] if p == args.process))
    nom = {ch: tpl.get(args.process, ch)[0] for ch in channels}
    nom_all = np.concatenate([nom[ch] for ch in channels])

    # other theory uncertainties: per bin, max(|up-1|, |down-1|) of each shape nuisance, in quadrature
    shapes = {n: e for n, (t, e) in card_nuisances(card, args.process, ["theo"]).items() if t == "shape"}
    band2 = np.zeros_like(nom_all)
    for name, entries in shapes.items():
        dev = []
        for ud in ("Up", "Down"):
            h = np.concatenate([tpl.get(args.process, ch, name + ud)[0] if ch in entries else nom[ch]
                                for ch in channels])
            scale = 1.0 if args.keep_rate else nom_all.sum() / h.sum()
            dev.append(np.abs(np.divide(h * scale, nom_all, out=np.ones_like(nom_all), where=nom_all > 0) - 1))
        band2 += np.maximum(*dev) ** 2
    band = np.sqrt(band2)

    with open(args.hdamp) as f:
        hd = json.load(f)
    n_nn = hd.get("n_nn") or {}
    cols = {}
    for v in ("nominal", "up", "down"):
        hs, es = [], []
        for ch in channels:
            h, e = hd["channels"][ch][v]
            h, e2 = to_card_bins(h, np.square(e), n_nn.get(ch, 6), len(nom[ch]))
            hs.append(h)
            es.append(np.sqrt(e2))
        cols[v] = (np.concatenate(hs), np.concatenate(es))
    hn, en = cols["nominal"]
    ratios = {}
    for v in ("up", "down"):
        hv, ev = cols[v]
        scale = 1.0 if args.keep_rate else hn.sum() / hv.sum()
        r = np.divide(hv * scale, hn, out=np.full_like(hn, np.nan), where=hn > 0)
        er = np.sqrt(np.divide(ev * scale, hn, out=np.zeros_like(hn), where=hn > 0) ** 2
                     + (r * np.divide(en, hn, out=np.zeros_like(hn), where=hn > 0)) ** 2)
        ratios[v] = (r, er, hv.sum() / hn.sum())
        inside = np.abs(r - 1) <= np.sqrt(band ** 2 + er ** 2)
        print(f"hdamp {v}: rate {hv.sum() / hn.sum():.3f}; {inside.sum()}/{len(r)} bins within the band "
              f"(band and hdamp stat. combined)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(channels), 1, figsize=(10, 3.0 * len(channels)), squeeze=False)
    start = 0
    for ax, ch in zip(axes[:, 0], channels):
        n = len(nom[ch])
        sl = slice(start, start + n)
        x = np.arange(n)
        edges = np.arange(n + 1) - 0.5
        b = band[sl]
        ax.fill_between(edges, np.append(1 - b, 1 - b[-1]), np.append(1 + b, 1 + b[-1]), step="post",
                        color="#bdbdbd", alpha=0.7, lw=0, label="other tt+B theory unc. (card, quadrature)")
        for k, (v, c, m) in enumerate((("up", "#1f4e9c", "o"), ("down", "#d9730d", "s"))):
            r, er, rate = ratios[v]
            ax.errorbar(x + (k - 0.5) * 0.25, r[sl], yerr=er[sl], fmt=m, ms=3.5, color=c, lw=1,
                        mfc=c if v == "up" else "white", label=f"hdamp {v} / nominal (sample stat.)")
        ax.axhline(1, color="black", lw=0.8)
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(0.5, 1.5)
        ax.set_ylabel("variation / nominal")
        ax.set_title(f"{args.process}, {ch}" + ("" if args.keep_rate else " (each variation scaled to the nominal total)"),
                     fontsize=9)
        ax.grid(alpha=0.3)
        start += n
    axes[0, 0].legend(fontsize=7, ncol=3, loc="upper left", frameon=False)
    axes[-1, 0].set_xlabel("datacard bin (NN bin x mass bin)")
    fig.tight_layout()
    stem = f"{args.outdir}/{args.process}_hdamp_vs_theory{'_rate' if args.keep_rate else ''}"
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)
    print(f"Wrote {stem}.png/.pdf")


if __name__ == "__main__":
    main()
