#!/usr/bin/env python3
"""
4FS ttbb cross sections normalized to the inclusive 5FS tt+B (Run 3 port of SemiLeptonic_2024/calc_ttbbXS.py,
the script the 2024 values of outputToPkl.XSEC_DICT came from, run there with a config without cuts).
Per decay (TTbb_X paired with TTToX):

    sigma_4FS = sigma_5FS * [N_5FS(tt+B) / N_5FS(all)] * [N_4FS(all) / N_4FS(tt+B)]

over all generated events, tt+B = genTtbarId % 100 >= 51, N = sum of sign(genWeight). The sums are the
theory_sums the processor fills before the skim (theory_norms.py: regions 'all' and 'ttB', variation
'nominal'), so the nominal outputs suffice and no special run is needed. The 4FS sums carry ttbb_scale_ref
(analysis nominal (mu0, mu0), custom_scale_factors), as the chain does (outputToPkl.get_norm).

    python3 calc_ttbb_xsec.py -i outputs/nominal/bkg

-i takes .coffea files or output dirs (merged/<dataset>.coffea or output_all.coffea).
"""
import argparse
import glob
import os

from analysis_years import find_era
from outputToPkl import XSEC_DICT, xsec_key

PAIRS = {"TTbb_SemiLeptonic": "TTToSemiLeptonic", "TTbb_Hadronic": "TTToHadronic", "TTbb_2L2Nu": "TTTo2L2Nu"}


def coffea_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            merged = sorted(glob.glob(os.path.join(p, "merged", "*.coffea")))
            files += merged or sorted(glob.glob(os.path.join(p, "output_all.coffea")))
        else:
            files.append(p)
    return files


def main():
    from coffea.util import load
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True)
    args = ap.parse_args()

    theory = {}
    for path in coffea_files(args.input):
        out = load(path)
        for sample, by_ds in out.get("theory_sums", {}).items():
            theory.setdefault(sample, {}).update(by_ds)
        del out
    if not theory:
        raise SystemExit("No theory_sums in the inputs (outputs made before the Run 3 datacard port?)")

    for s4, s5 in PAIRS.items():
        for d4, t4 in theory.get(s4, {}).items():
            d5s = [d for d in theory.get(s5, {}) if find_era(d) == find_era(d4)]
            if not d5s:
                print(f"  {d4}: no {s5} theory_sums for its era, skipped")
                continue
            t5 = theory[s5][d5s[0]]
            f5 = t5["ttB"]["nominal"] / t5["all"]["nominal"]
            f4 = t4["ttB"]["nominal"] / t4["all"]["nominal"]
            sigma5 = XSEC_DICT[xsec_key(s5, d5s[0])]
            sigma4 = sigma5 * f5 / f4
            used = XSEC_DICT[xsec_key(s4, d4)]
            print(f"  {d4:22s} {sigma4:9.4f} pb   (XSEC_DICT {used:9.4f} pb, {sigma4 / used - 1:+.1%})   "
                  f"tt+B fraction 5FS {f5:.5f} (sigma_5FS {sigma5} pb), 4FS {f4:.4f}")


if __name__ == "__main__":
    main()
