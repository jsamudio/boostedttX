#!/usr/bin/env python3
"""
GloParT XbbVsQCD working points of the Run 3 analysis from the mistag rate of light-flavour fatjets
in tt+LF events (output of the xbbwp group):

    python3 run_pipeline.py inference --backend pbs --groups xbbwp
    python3 derive_xbb_wps.py [-i outputs/nominal/xbbwp/output_all.coffea] [--plots plots/xbbwp]

Mistag fatjets: all FatJetGood of the analysis preselection (1 lepton, >= 5 AK4, >= 1 AK8; add
--min-bjets 2 for the btag_mask category) in tt+LF events (TTToSemiLeptonic, TTTo2L2Nu) with no B or
C hadron (nBHadrons == nCHadrons == 0). The cut of each WP is where the weighted fraction of these
fatjets above it equals the target rate. Signal efficiency: fatjets containing both b quarks of the
H/Z (n_boson_b == 2) in ttH(bb) / ttZ(bb). The result goes to analysis_years.XBB_WP_DICT.
"""
import argparse
import os

import numpy as np
from coffea.util import load

HERE = os.path.dirname(os.path.abspath(__file__))
CATEGORY = "presel"
TARGETS = {"L": 0.10, "M": 0.05, "T": 0.01}
MISTAG_KEYS = ["TTToSemiLeptonic__tt+LF", "TTTo2L2Nu__tt+LF"]
SIGNAL_KEYS = ["ttHTobb", "TTZToQQ__Zbb"]
PT_BINS = [200, 300, 400, 600, np.inf]


def fatjets(output, keys, min_bjets=0):
    """Per-fatjet arrays (concatenated over keys / datasets) with the event weight normalized by the
    sum of genweights, plus the per-event leading Xbb (the ZH candidate)."""
    out = {k: [] for k in ("pt", "mass", "xbb", "nB", "nC", "n_boson_b", "weight")}
    cand = {"xbb": [], "weight": [], "light": []}
    for key in keys:
        if key not in output["columns"]:
            print(f"  [skip] {key} not in the output")
            continue
        for dataset, cats in output["columns"][key].items():
            c = {name: col.value for name, col in cats[CATEGORY].items()}
            n = c["FatJetGood_N"]
            sel_evt = c["events_nbJetGood"] >= min_bjets
            w_evt = c["weight"] / output["sum_genweights"][dataset]
            sel = np.repeat(sel_evt, n)
            w = np.repeat(w_evt, n)[sel]
            out["pt"].append(c["FatJetGood_pt"][sel])
            out["mass"].append(c["FatJetGood_mass"][sel])
            out["xbb"].append(c["FatJetGood_btagBB"][sel])
            out["nB"].append(c["FatJetGood_nBHadrons"][sel])
            out["nC"].append(c["FatJetGood_nCHadrons"][sel])
            out["n_boson_b"].append(c["FatJetGood_n_boson_b"][sel])
            out["weight"].append(w)
            # candidate = leading Xbb fatjet of each event
            offsets = np.concatenate([[0], np.cumsum(n)])
            xbb_all = c["FatJetGood_btagBB"]
            light_all = (c["FatJetGood_nBHadrons"] == 0) & (c["FatJetGood_nCHadrons"] == 0)
            for i in np.nonzero(sel_evt & (n > 0))[0]:
                j = offsets[i] + np.argmax(xbb_all[offsets[i]:offsets[i + 1]])
                cand["xbb"].append(xbb_all[j])
                cand["light"].append(light_all[j])
                cand["weight"].append(w_evt[i])
    out = {k: np.concatenate(v) if v else np.array([]) for k, v in out.items()}
    return out, {k: np.asarray(v) for k, v in cand.items()}


def cut_at_rate(score, weight, rate):
    """Smallest cut with weighted fraction(score > cut) <= rate."""
    order = np.argsort(score)[::-1]
    frac = np.cumsum(weight[order]) / weight.sum()
    idx = np.searchsorted(frac, rate, side="right")
    return float(score[order][min(idx, len(order) - 1)])


def frac_above(score, weight, cut):
    return float(weight[score > cut].sum() / weight.sum()) if len(score) else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", default=os.path.join(HERE, "outputs/nominal/xbbwp/output_all.coffea"))
    ap.add_argument("--era", default="2024")
    ap.add_argument("--min-bjets", type=int, default=0, help="require nbJetGood >= N (2 = btag_mask)")
    ap.add_argument("--all-flavours", action="store_true",
                    help="mistag on all tt+LF fatjets, not only those without B/C hadrons")
    ap.add_argument("--plots", default=None, help="directory for the Xbb distributions / ROC")
    args = ap.parse_args()

    output = load(args.input)
    bkg, bkg_cand = fatjets(output, MISTAG_KEYS, args.min_bjets)
    sig, _ = fatjets(output, SIGNAL_KEYS, args.min_bjets)
    light = np.ones(len(bkg["xbb"]), dtype=bool) if args.all_flavours else (bkg["nB"] == 0) & (bkg["nC"] == 0)
    bb = sig["n_boson_b"] == 2
    print(f"mistag fatjets: {light.sum()} ({'all flavours' if args.all_flavours else 'light'}) of {len(light)}; "
          f"signal H/Z->bb fatjets: {bb.sum()} of {len(bb)}")
    if light.sum() == 0:
        raise SystemExit("No mistag fatjets found")

    xb, wb = bkg["xbb"][light], bkg["weight"][light]
    wps = {wp: cut_at_rate(xb, wb, rate) for wp, rate in TARGETS.items()}

    print(f"\n{'WP':>3} {'target':>7} {'cut':>8} {'mistag':>8} {'sig eff':>8} {'cand. mistag':>13}")
    cand_light = bkg_cand["light"] | args.all_flavours
    for wp, cut in wps.items():
        print(f"{wp:>3} {TARGETS[wp]:>7.3f} {cut:>8.4f} {frac_above(xb, wb, cut):>8.4f} "
              f"{frac_above(sig['xbb'][bb], sig['weight'][bb], cut):>8.3f} "
              f"{frac_above(bkg_cand['xbb'][cand_light], bkg_cand['weight'][cand_light], cut):>13.4f}")
    print("  (cand. mistag: leading-Xbb fatjet of the event, light flavour)")

    print("\nCut per fatjet pT bin (cross-check of the pT dependence):")
    print(f"{'pT bin':>12} " + " ".join(f"{wp:>8}" for wp in TARGETS) + f" {'n jets':>8}")
    pt = bkg["pt"][light]
    for lo, hi in zip(PT_BINS[:-1], PT_BINS[1:]):
        m = (pt >= lo) & (pt < hi)
        if m.sum() < 20:
            continue
        cuts = [cut_at_rate(xb[m], wb[m], r) for r in TARGETS.values()]
        print(f"{f'{lo}-{hi}':>12} " + " ".join(f"{c:>8.4f}" for c in cuts) + f" {m.sum():>8}")

    print("\nanalysis_years.XBB_WP_DICT entry:")
    keys = {"L": "LP", "M": "MP", "T": "HP"}          # XBB_WP_DICT naming (as in Run 2)
    print(f'    "{args.era}": {{' + ", ".join(f'"{keys[wp]}": {cut:.4f}' for wp, cut in wps.items()) + "},")

    if args.plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(args.plots, exist_ok=True)
        bins = np.linspace(0, 1, 51)
        fig, ax = plt.subplots()
        ax.hist(xb, bins, weights=wb / wb.sum(), histtype="step", label="tt+LF light fatjets")
        if bb.any():
            ax.hist(sig["xbb"][bb], bins, weights=sig["weight"][bb] / sig["weight"][bb].sum(),
                    histtype="step", label="H/Z->bb fatjets")
        for wp, cut in wps.items():
            ax.axvline(cut, ls="--", c="gray")
            ax.text(cut, ax.get_ylim()[1] * 0.9, wp)
        ax.set_yscale("log"); ax.set_xlabel("GloParT XbbVsQCD"); ax.set_ylabel("fraction"); ax.legend()
        fig.savefig(os.path.join(args.plots, f"xbb_{args.era}.pdf"))
        if bb.any():
            cuts = np.linspace(0, 1, 201)
            fig, ax = plt.subplots()
            ax.plot([frac_above(sig["xbb"][bb], sig["weight"][bb], c) for c in cuts],
                    [frac_above(xb, wb, c) for c in cuts])
            for wp, cut in wps.items():
                ax.plot(frac_above(sig["xbb"][bb], sig["weight"][bb], cut), frac_above(xb, wb, cut), "o", label=wp)
            ax.set_yscale("log"); ax.set_xlabel("H/Z->bb efficiency"); ax.set_ylabel("light mistag rate")
            ax.legend()
            fig.savefig(os.path.join(args.plots, f"xbb_roc_{args.era}.pdf"))
        print(f"plots in {args.plots}")


if __name__ == "__main__":
    main()
