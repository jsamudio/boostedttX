"""
Generator-level bb pair of the EFT tt+bb sample, and the EFT-coefficient sums of its pT spectrum
over ALL generated events (before the skim), for plot_eft_pt.py.

gen bb pair: the two hard-process b quarks that come from neither a top nor a W/Z/H, i.e. the extra
g -> bb (provenance 1 of ../common/parton_matching.hard_partons with extra_b_as_boson=True, at the
first copy ~ matrix-element level). gen_bb_pair() returns pT, mass and dR of the pair (-1 if the
event does not have exactly two).

Accumulator (fill_eft_genbb, called from spanet_inference.process_extra_before_skim for EFT samples):

    output["eft_genbb"][dataset][meta] = array (len(GEN_BB_PT_BINS) - 1, n_coefficients)

the per-bin sum of the EFTfitCoefficients vectors, gen bb pT clipped into the last bin; events
without a gen bb pair are not counted. meta is a JSON string with the coefficient structure
(index1, index2, WC names), so that the chunks are summed per structure by the accumulator merge.
For a WC point c (c_0 = 1 for the SM): spectrum(c) = sums @ (c[index1] * c[index2]).
"""
import json

import awkward as ak
import numpy as np


GEN_BB_PT_BINS = np.arange(0, 1025, 25)
EFT_COLUMNS = ("EFTfitCoefficients", "EFTfitCoefficientIndex1", "EFTfitCoefficientIndex2", "WCnames")


def gen_bb_pair(events):
    from parton_matching import hard_partons      # configs/common (on the path in the processor)
    partons = hard_partons(events.GenPart, "isFirstCopy", extra_b_as_boson=True)
    b = partons[(abs(partons.pdgId) == 5) & (partons.provenance == 1)]
    two = ak.to_numpy(ak.num(b) == 2)
    b = ak.pad_none(b, 2, clip=True)
    pt1, pt2 = ak.to_numpy(ak.fill_none(b.pt[:, 0], 0.)), ak.to_numpy(ak.fill_none(b.pt[:, 1], 0.))
    eta1, eta2 = ak.to_numpy(ak.fill_none(b.eta[:, 0], 0.)), ak.to_numpy(ak.fill_none(b.eta[:, 1], 0.))
    phi1, phi2 = ak.to_numpy(ak.fill_none(b.phi[:, 0], 0.)), ak.to_numpy(ak.fill_none(b.phi[:, 1], 0.))
    m1, m2 = ak.to_numpy(ak.fill_none(b.mass[:, 0], 0.)), ak.to_numpy(ak.fill_none(b.mass[:, 1], 0.))
    px = pt1 * np.cos(phi1) + pt2 * np.cos(phi2)
    py = pt1 * np.sin(phi1) + pt2 * np.sin(phi2)
    pz = pt1 * np.sinh(eta1) + pt2 * np.sinh(eta2)
    energy = (np.sqrt((pt1 * np.cosh(eta1)) ** 2 + m1 ** 2) + np.sqrt((pt2 * np.cosh(eta2)) ** 2 + m2 ** 2))
    pt = np.hypot(px, py)
    mass = np.sqrt(np.maximum(energy ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.))
    dphi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
    dr = np.hypot(eta1 - eta2, dphi)
    return tuple(np.where(two, x, -1.) for x in (pt, mass, dr))


def is_eft_sample(events):
    return all(c in events.fields for c in EFT_COLUMNS)


def coefficient_structure(events):
    """JSON string of (index1, index2, WC names as stored) of the first event."""
    return json.dumps({
        "index1": [int(i) for i in ak.to_numpy(events.EFTfitCoefficientIndex1[0])],
        "index2": [int(i) for i in ak.to_numpy(events.EFTfitCoefficientIndex2[0])],
        "wcnames": [int(v) for v in ak.to_numpy(events.WCnames[0])],
    })


def fill_eft_genbb(events, output, dataset):
    """Add this chunk's coefficient sums per gen bb pT bin to output['eft_genbb'][dataset]."""
    if len(events) == 0 or not is_eft_sample(events):
        return
    pt, _, _ = gen_bb_pair(events)
    coefs = ak.to_numpy(ak.to_regular(events.EFTfitCoefficients, axis=1)).astype(np.float64)
    sums = np.zeros((len(GEN_BB_PT_BINS) - 1, coefs.shape[1]))
    ok = pt >= 0
    idx = np.clip(np.digitize(pt[ok], GEN_BB_PT_BINS) - 1, 0, len(GEN_BB_PT_BINS) - 2)
    np.add.at(sums, idx, coefs[ok])
    per_dataset = output.setdefault("eft_genbb", {}).setdefault(dataset, {})
    key = coefficient_structure(events)
    per_dataset[key] = per_dataset[key] + sums if key in per_dataset else sums
    totals = output.setdefault("eft_genbb_total", {}).setdefault(dataset, {})
    total = coefs.sum(axis=0)
    totals[key] = totals[key] + total if key in totals else total
