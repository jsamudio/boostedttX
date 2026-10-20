"""
Sums of the theory weights over ALL generated events (before the skim), used by make_process_norms.py
to make the theory up/down templates keep their acceptance effect but not the change of the total
cross section (the method of the previous analysis; see make_process_norms.py).

Called once per chunk from the processor (spanet_inference.py, process_extra_before_skim):

    output["theory_sums"][sample][dataset][region][variation] = sum of sign(genWeight) * ratio

  variation : "nominal" (ratio = 1) and the theory weight columns used by makeDatacard.py:
              mu_r_up/down, mu_f_up/down, mu_rf_up/down, isr_up/down, fsr_up/down,
              pdf_0 .. pdf_<n-1> (one per LHEPdfWeight member; n identifies the PDF set, see
              pdf_uncertainty.py)
  region    : "all"      every event of the dataset
              "ttB"      ttbar samples (5FS TTTo, 4FS TTbb): events with genTtbarId % 100 >= 51 (tt+B)
              "genbinN"  signal samples: generator-level Z/H pT bin N (bins of GEN_PT_BINS) with
                         |y(Z/H)| < SIGNAL_MAX_ABS_Y; TTZToQQ events with Z -> bb are left out
                         (they are not part of the ttZ process, see genmatcher.py: 'old_ttZbb')

The ratios are taken with the SAME functions as the event weights in the processor
(sf_scaleweights, sf_partonshower_isr/fsr, sf_pdfweights), so the sums match the templates exactly.
"""
import awkward as ak
import numpy as np

from custom_scale_factors import sf_scaleweights, sf_partonshower_isr, sf_partonshower_fsr, sf_pdfweights

# Generator-level Z/H pT bins of the signal processes (ttH0..ttH3, ttZ0..ttZ3), in GeV
GEN_PT_BINS = [0, 200, 300, 450, np.inf]
# STXS rapidity requirement on the generator-level Z/H; None = no requirement
SIGNAL_MAX_ABS_Y = 2.5

SIGNAL_SAMPLES = ["ttHTobb", "ttHToNonbb", "TTZToBB", "TTZToQQ", "TTZToLLNuNu"]
TTBAR_SAMPLES_PREFIX = ("TTTo", "TTbb")


def theory_ratios(events):
    """{variation: per-event ratio to the nominal weight}; 'nominal' is 1."""
    n = len(events)
    ratios = {"nominal": np.ones(n)}
    scale = sf_scaleweights(events)          # None if the sample has no LHEScaleWeight
    if scale is not None:
        names = ["mu_r_up", "mu_r_down", "mu_f_up", "mu_f_down", "mu_rf_up", "mu_rf_down"]
        for name, values in zip(names, scale):
            ratios[name] = ak.to_numpy(values)
    if "PSWeight" in events.fields and ak.min(ak.num(events.PSWeight)) >= 4:
        _, ratios["isr_up"], ratios["isr_down"] = (ak.to_numpy(x) for x in sf_partonshower_isr(events))
        _, ratios["fsr_up"], ratios["fsr_down"] = (ak.to_numpy(x) for x in sf_partonshower_fsr(events))
    for i, values in enumerate(sf_pdfweights(events)):
        ratios[f"pdf_{i}"] = values
    return ratios


def gen_boson(events):
    """pT, rapidity and Z->bb flag of the generator-level Z/H (status 62, as in genmatcher.py).
    Events without such a boson get pT = -1."""
    gen = events.GenPart
    is_boson = ((gen.pdgId == 23) | (gen.pdgId == 25)) & (gen.status == 62)
    boson = ak.firsts(gen[is_boson])
    pt = ak.to_numpy(ak.fill_none(boson.pt, -1.0))
    eta = ak.to_numpy(ak.fill_none(boson.eta, 0.0))
    mass = ak.to_numpy(ak.fill_none(boson.mass, 0.0))
    # rapidity from pT, eta and mass: y = 0.5 * ln((E + pz) / (E - pz))
    pz = pt * np.sinh(eta)
    energy = np.sqrt((pt * np.cosh(eta)) ** 2 + mass ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        rapidity = np.where(pt > 0, 0.5 * np.log((energy + pz) / (energy - pz)), 99.0)
    # Z -> bb: two b quarks whose mother is a status-62 Z (same definition as genmatcher.py)
    mother_id = gen.pdgId[gen.genPartIdxMother]
    mother_status = gen.status[gen.genPartIdxMother]
    b_from_z = (abs(gen.pdgId) == 5) & (mother_id == 23) & (mother_status == 62) & (gen.genPartIdxMother >= 0)
    is_zbb = ak.to_numpy(ak.sum(b_from_z, axis=1) == 2)
    return pt, rapidity, is_zbb


def region_masks(events, sample):
    """{region: boolean mask} for this sample (see the module docstring)."""
    n = len(events)
    masks = {"all": np.ones(n, dtype=bool)}
    if sample.startswith(TTBAR_SAMPLES_PREFIX):
        masks["ttB"] = ak.to_numpy(events.genTtbarId % 100 >= 51)
    if sample in SIGNAL_SAMPLES:
        pt, rapidity, is_zbb = gen_boson(events)
        keep = pt >= 0
        if SIGNAL_MAX_ABS_Y is not None:
            keep = keep & (np.abs(rapidity) < SIGNAL_MAX_ABS_Y)
        if sample == "TTZToQQ":
            keep = keep & ~is_zbb
        for i in range(len(GEN_PT_BINS) - 1):
            in_bin = (pt >= GEN_PT_BINS[i]) & (pt < GEN_PT_BINS[i + 1])
            masks[f"genbin{i}"] = keep & in_bin
    return masks


def fill_theory_sums(events, output, sample, dataset):
    """Add this chunk's sums to output['theory_sums'][sample][dataset][region][variation]."""
    if len(events) == 0:
        return
    sign = np.sign(ak.to_numpy(events.genWeight))
    ratios = theory_ratios(events)
    out = output.setdefault("theory_sums", {}).setdefault(sample, {}).setdefault(dataset, {})
    for region, mask in region_masks(events, sample).items():
        sums = out.setdefault(region, {})
        for variation, ratio in ratios.items():
            sums[variation] = sums.get(variation, 0.0) + float(np.sum(sign[mask] * ratio[mask]))
