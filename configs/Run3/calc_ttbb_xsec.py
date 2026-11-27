#!/usr/bin/env python3
"""
4FS ttbb cross sections normalized to the 5FS tt+B (Run 3 port of SemiLeptonic_2024/calc_ttbbXS.py, the
script the 2024 values of outputToPkl.XSEC_DICT came from). Per decay (TTbb_X paired with TTToX):

    sigma_4FS = sigma_5FS * [N_5FS(tt+B, selected) / N_5FS(all generated)]
                          * [N_4FS(all generated) / N_4FS(tt+B, selected)]

so that the 4FS tt+B yield after the selection equals the 5FS tt+B yield. As the old script, N = sum of
sign(genWeight); the 4FS sums carry ttbb_scale_ref (analysis nominal (mu0, mu0), custom_scale_factors), as the
chain does (the all-generated 4FS sum is theory_sums[...]['all']['nominal'], the denominator of
outputToPkl.get_norm). 'Selected' = the tt+B subsample in the btag_mask category of the nominal output
(the old script: its 'baseline' category); with --datacard-cuts also makeDatacard.cuts (Xbb >= LP, ...).
The inclusive version (all generated events, genTtbarId % 100 >= 51, from the theory_sums) is printed too.

    python3 calc_ttbb_xsec.py -i outputs/nominal/bkg [--datacard-cuts]

-i takes .coffea files or output dirs (merged/<dataset>.coffea or output_all.coffea).
"""
import argparse
import glob
import os

import numpy as np

from analysis_years import find_era, bbvl_wp
from outputToPkl import XSEC_DICT, xsec_key

CATEGORY = "btag_mask"
PAIRS = {"TTbb_SemiLeptonic": "TTToSemiLeptonic", "TTbb_Hadronic": "TTToHadronic", "TTbb_2L2Nu": "TTTo2L2Nu"}
CUT_VARS = ["n_ak4jets", "n_b_outZH", "ZH_bbvLscore", "MET_pt", "ZH_M"]


def coffea_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            merged = sorted(glob.glob(os.path.join(p, "merged", "*.coffea")))
            files += merged or sorted(glob.glob(os.path.join(p, "output_all.coffea")))
        else:
            files.append(p)
    return files


def selected_sum(cols, era, with_ref, datacard_cuts):
    """sum of sign(genWeight) [x ttbb_scale_ref] over the events of one dataset's columns."""
    w = np.sign(cols["events_genWeight"].value)
    if with_ref:
        w = w * cols["events_ttbb_scale_ref"].value
    if datacard_cuts:
        v = {k: cols[f"events_{k}"].value for k in CUT_VARS}
        mask = ((v["n_ak4jets"] >= 5) & (v["n_b_outZH"] == 2) & (v["ZH_bbvLscore"] >= bbvl_wp(era))
                & (v["MET_pt"] >= 20) & (v["ZH_M"] >= 50) & (v["ZH_M"] <= 200))
        w = w[mask]
    return float(np.sum(w))


def main():
    from coffea.util import load
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True)
    ap.add_argument("--datacard-cuts", action="store_true", help="also apply makeDatacard.cuts to the selected events")
    args = ap.parse_args()

    columns, genw, theory = {}, {}, {}
    for path in coffea_files(args.input):
        out = load(path)
        for sub, by_ds in out.get("columns", {}).items():
            if sub.endswith("__tt+B"):
                for ds, cats in by_ds.items():
                    if CATEGORY in cats:
                        columns.setdefault(sub, {})[ds] = cats[CATEGORY]["nominal"]
        genw.update(out.get("sum_signOf_genweights", {}))
        for sample, by_ds in out.get("theory_sums", {}).items():
            theory.setdefault(sample, {}).update(by_ds)
        del out

    print(f"Selection: tt+B subsample, {CATEGORY}" + (" + makeDatacard.cuts" if args.datacard_cuts else ""))
    for s4, s5 in PAIRS.items():
        for d4, cols4 in columns.get(f"{s4}__tt+B", {}).items():
            era = find_era(d4)
            d5s = [d for d in columns.get(f"{s5}__tt+B", {}) if find_era(d) == era]
            if not d5s:
                print(f"  {d4}: no {s5} tt+B columns for {era}, skipped")
                continue
            d5 = d5s[0]
            sigma5 = XSEC_DICT[xsec_key(s5, d5)]
            n5_sel = selected_sum(columns[f"{s5}__tt+B"][d5], era, False, args.datacard_cuts)
            n5_all = genw[d5]
            n4_sel = selected_sum(cols4, era, True, args.datacard_cuts)
            try:
                n4_all = theory[s4][d4]["all"]["nominal"]
            except KeyError:
                raise SystemExit(f"No theory_sums for {s4}/{d4}: the outputs predate the Run 3 datacard port")
            sigma4 = sigma5 * (n5_sel / n5_all) * (n4_all / n4_sel)
            used = XSEC_DICT[xsec_key(s4, d4)]
            line = (f"  {d4:22s} selected match {sigma4:9.4f} pb   XSEC_DICT {used:9.4f} pb ({sigma4 / used - 1:+.1%})"
                    f"   [5FS tt+B {n5_sel / n5_all:.5f} of all; 4FS tt+B {n4_sel / n4_all:.4f} of all]")
            t5, t4 = theory.get(s5, {}).get(d5), theory[s4][d4]
            if t5 and "ttB" in t5 and "ttB" in t4:
                incl = (sigma5 * (t5["ttB"]["nominal"] / t5["all"]["nominal"])
                        / (t4["ttB"]["nominal"] / t4["all"]["nominal"]))
                line += f"\n  {'':22s} inclusive match {incl:8.4f} pb (genTtbarId % 100 >= 51, all generated)"
            print(line)


if __name__ == "__main__":
    main()
