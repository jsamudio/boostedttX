#!/usr/bin/env python3
"""
Write process_norms/process_norms_<year>.json: the factors makeDatacard.py multiplies the theory
up/down templates by, so that a variation keeps its acceptance/shape effect but not the change of
the total cross section (that part is covered by the lnN uncertainties, and for tt+B by its
floating rate). Method of the previous analysis.

Input: the theory-weight sums over ALL generated events, accumulated by the processor
(theory_norms.py) in output['theory_sums'] of the nominal .coffea outputs.

    python3 make_process_norms.py -i outputs/nominal/signal/output_all.coffea \\
                                     outputs/nominal/bkg/output_all.coffea

For a variation v (e.g. mu_r_up) and N_x = sum of sign(genWeight) * ratio_x over the events:

  ttH0..ttH3, ttZ0..ttZ3  N_nom / N_v  in the generator-level Z/H pT bin (|y| < 2.5)
                          (mu_r, mu_f, mu_rf, isr, fsr, pdf, alphas)
  TTBar                   N_nom / N_v  over the full 5FS sample (mu_r, mu_f, mu_rf, pdf, alphas;
                          ISR/FSR of TTBar are not normalized, see makeDatacard.py)
  ttbb (= tt_B)           paired 4FS (TTbb_X) and 5FS (TTToX) samples:
      mu_r, mu_f, mu_rf : (4FS ttB_nom / 4FS ttB_v) * (5FS ttB_v / 5FS all_v) / (5FS ttB_nom / 5FS all_nom)
                          = remove the 4FS change, apply the change of the tt+B FRACTION in the 5FS
      isr, fsr          : (5FS ttB_v / 5FS ttB_nom) / (4FS ttB_v / 4FS ttB_nom)
                          = remove the 4FS change, apply the change of the tt+B YIELD in the 5FS
      pdf               : N_nom / N_v of the 4FS tt+B only (4FS PDF set, no alpha_s members)

  pdf / alphas up and down are computed from the sums of the LHEPdfWeight members with the same
  functions as the templates (pdf_uncertainty.py), so the factors remove exactly their rate part.

ttX (no theory reference cross section for the mix of samples): instead of factors, the
uncertainty of the total ttX cross section from the samples' own weights, used by makeDatacard.py
as lnN values (previous analysis: "envelope of renorm./fact. scales", "norm. uncertainty due to PDF"):
  ttX_xsec_unc  scale_up / scale_down  max / min over mu_r, mu_f, mu_rf up/down of N_v / N_nom
                pdf_alphas             1 + (relative PDF (+) alpha_s uncertainty, in quadrature)
  (per variation, the xsec x lumi weighted average over the ttX datasets = the change of the total)

A process made of several datasets (e.g. ttZ = TTZToBB + TTZToQQ + TTZToLLNuNu, or the 2016 eras)
gets the average of the dataset factors, weighted by cross section x luminosity (as before).
"""
import argparse
import json
import os

import numpy as np

from coffea.util import load

from analysis_years import find_era, ERA_TO_YEAR, LUMI_FB
from outputToPkl import XSEC_DICT, ERA_SUFFIX
from theory_norms import GEN_PT_BINS
from pdf_uncertainty import pdf_up_down, alphas_up_down, has_alphas

SCALE_VARIATIONS = ["mu_r_up", "mu_r_down", "mu_f_up", "mu_f_down", "mu_rf_up", "mu_rf_down"]
PS_VARIATIONS = ["isr_up", "isr_down", "fsr_up", "fsr_down"]

TTX_SAMPLE = "ttX"

SIGNAL_PROCESSES = {
    "ttH": ["ttHTobb", "ttHToNonbb"],
    "ttZ": ["TTZToBB", "TTZToQQ", "TTZToLLNuNu"],
}
TTBAR_5FS = ["TTToSemiLeptonic", "TTTo2L2Nu", "TTToHadronic"]
TTBB_4FS_TO_5FS = {                          # each 4FS sample with the 5FS sample of the same decay
    "TTbb_SemiLeptonic": "TTToSemiLeptonic",
    "TTbb_2L2Nu": "TTTo2L2Nu",
    "TTbb_Hadronic": "TTToHadronic",
}
# TTZToQQ events with Z -> bb are not part of ttZ (TTZToBB covers them): use the remaining cross section
XSEC_OVERRIDE = {"TTZToQQ": XSEC_DICT["TTZToQQ"] - XSEC_DICT["TTZToBB"]}

OUT_DIR = "process_norms"


# ---------------------------------------------------------------- helpers
def dataset_weight(sample, dataset):
    """Share of a dataset in its process: cross section x luminosity of its era."""
    name = ERA_SUFFIX.sub("", dataset)
    xsec = XSEC_OVERRIDE.get(sample, XSEC_DICT.get(name, XSEC_DICT.get(sample)))
    if xsec is None:
        raise KeyError(f"No cross section for {dataset} (sample {sample}) in outputToPkl.XSEC_DICT")
    return xsec * LUMI_FB[find_era(dataset)]


def year_of(dataset):
    return ERA_TO_YEAR[find_era(dataset)]


def ratio(sums, region, variation):
    """N_nom / N_variation in one region of one dataset."""
    if variation not in sums[region]:
        raise KeyError(f"Variation '{variation}' missing in region '{region}'")
    if sums[region][variation] == 0:
        raise ZeroDivisionError(f"No events in region '{region}' for variation '{variation}'")
    return sums[region]["nominal"] / sums[region][variation]


def weighted_average(entries):
    """entries: list of (weight, {variation: factor}) -> {variation: weighted average}."""
    keys = set(entries[0][1])
    if any(set(f) != keys for _, f in entries[1:]):
        raise KeyError("Datasets of one process have different variations "
                       "(e.g. PDF weights in some datasets only): " + str([sorted(f) for _, f in entries]))
    total = sum(w for w, _ in entries)
    return {v: sum(w * f[v] for w, f in entries) / total for v in entries[0][1]}


def pdf_members(sums, region):
    """Sums of the LHEPdfWeight members in a region, in member order (empty without PDF weights)."""
    n = sum(1 for v in sums[region] if v.startswith("pdf_"))
    return np.array([sums[region][f"pdf_{i}"] for i in range(n)])


def pdf_factors(sums, region):
    """N_nom / N_v for v = pdf_up, pdf_down (+ alphas_up, alphas_down for the Hessian set).
    Empty if the dataset has no PDF weights (outputs made before they were added)."""
    members = pdf_members(sums, region)
    if len(members) == 0:
        return {}
    nominal = sums[region]["nominal"]
    up, down = pdf_up_down(members)
    factors = {"pdf_up": nominal / up, "pdf_down": nominal / down}
    if has_alphas(len(members)):
        up, down = alphas_up_down(members)
        factors.update({"alphas_up": nominal / up, "alphas_down": nominal / down})
    return factors


# ---------------------------------------------------------------- per process type
def signal_factors(theory_sums, samples, year, gen_bin):
    entries = []
    for sample in samples:
        for dataset, sums in theory_sums.get(sample, {}).items():
            if year_of(dataset) != year:
                continue
            region = f"genbin{gen_bin}"
            factors = {v: ratio(sums, region, v) for v in SCALE_VARIATIONS + PS_VARIATIONS}
            factors.update(pdf_factors(sums, region))
            entries.append((dataset_weight(sample, dataset), factors))
    return weighted_average(entries) if entries else None


def ttbar_factors(theory_sums, year):
    entries = []
    for sample in TTBAR_5FS:
        for dataset, sums in theory_sums.get(sample, {}).items():
            if year_of(dataset) != year:
                continue
            factors = {v: ratio(sums, "all", v) for v in SCALE_VARIATIONS}
            factors.update(pdf_factors(sums, "all"))
            entries.append((dataset_weight(sample, dataset), factors))
    return weighted_average(entries) if entries else None


def ttbb_factors(theory_sums, year):
    entries = []
    for sample_4fs, sample_5fs in TTBB_4FS_TO_5FS.items():
        for dataset_4fs, s4 in theory_sums.get(sample_4fs, {}).items():
            if year_of(dataset_4fs) != year:
                continue
            era = find_era(dataset_4fs)
            matches = [d for d in theory_sums.get(sample_5fs, {}) if find_era(d) == era]
            if not matches:
                raise KeyError(f"No {sample_5fs} dataset for era {era} (needed for {dataset_4fs})")
            s5 = theory_sums[sample_5fs][matches[0]]
            factors = {}
            for v in SCALE_VARIATIONS:
                fraction_5fs_var = s5["ttB"][v] / s5["all"][v]
                fraction_5fs_nom = s5["ttB"]["nominal"] / s5["all"]["nominal"]
                factors[v] = ratio(s4, "ttB", v) * fraction_5fs_var / fraction_5fs_nom
            for v in PS_VARIATIONS:
                change_5fs = s5["ttB"][v] / s5["ttB"]["nominal"]
                change_4fs = s4["ttB"][v] / s4["ttB"]["nominal"]
                factors[v] = change_5fs / change_4fs
            factors.update(pdf_factors(s4, "ttB"))
            entries.append((dataset_weight(sample_4fs, dataset_4fs), factors))
    return weighted_average(entries) if entries else None


def ttx_xsec_uncertainty(theory_sums, year):
    """lnN values of the total ttX cross section (see the module docstring)."""
    entries = []
    for dataset, sums in theory_sums.get(TTX_SAMPLE, {}).items():
        if year_of(dataset) != year:
            continue
        changes = {v: 1 / ratio(sums, "all", v) for v in SCALE_VARIATIONS}      # N_v / N_nom
        members = pdf_members(sums, "all")
        if len(members):
            changes["pdf"] = pdf_up_down(members)[0] / members[0] - 1
            changes["alphas"] = alphas_up_down(members)[0] / members[0] - 1 if has_alphas(len(members)) else 0.0
        entries.append((dataset_weight(TTX_SAMPLE, dataset), changes))
    if not entries:
        return None
    average = weighted_average(entries)
    unc = {"scale_up": max(average[v] for v in SCALE_VARIATIONS),
           "scale_down": min(average[v] for v in SCALE_VARIATIONS)}
    if "pdf" in average:
        unc["pdf_alphas"] = 1 + float(np.hypot(average["pdf"], average["alphas"]))
    else:
        print(f"  [WARNING] {year}: no PDF weights in the ttX sums, ttX PDF uncertainty not derived")
    return unc


# ---------------------------------------------------------------- main
def read_theory_sums(paths):
    theory_sums = {}
    for path in paths:
        print(f"Reading {path}")
        out = load(path)
        if "theory_sums" not in out:
            print(f"  [WARNING] no 'theory_sums' in {path}: produced before theory_norms.py was added?")
        theory_sums.update(out.get("theory_sums", {}))
        del out
    return theory_sums


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True, help="nominal .coffea outputs (signal, bkg)")
    ap.add_argument("-o", "--outdir", default=OUT_DIR)
    args = ap.parse_args()

    theory_sums = read_theory_sums(args.input)
    years = sorted({year_of(d) for datasets in theory_sums.values() for d in datasets})
    os.makedirs(args.outdir, exist_ok=True)

    for year in years:
        norms = {}
        for process, samples in SIGNAL_PROCESSES.items():
            for gen_bin in range(len(GEN_PT_BINS) - 1):
                factors = signal_factors(theory_sums, samples, year, gen_bin)
                if factors:
                    norms[f"{process}{gen_bin}"] = factors
        for name, factors in [("TTBar", ttbar_factors(theory_sums, year)),
                              ("ttbb", ttbb_factors(theory_sums, year)),
                              ("ttX_xsec_unc", ttx_xsec_uncertainty(theory_sums, year))]:
            if factors:
                norms[name] = factors

        path = os.path.join(args.outdir, f"process_norms_{year}.json")
        with open(path, "w") as f:
            json.dump({year: norms}, f, indent=2)
        print(f"\nWrote {path}")
        for process, factors in norms.items():
            print(f"  {process:12} " + "  ".join(f"{v}={x:.3f}" for v, x in factors.items()))


if __name__ == "__main__":
    main()
