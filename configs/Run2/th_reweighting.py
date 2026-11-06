"""
SM reweighting of the tH samples (Run 2 UL THQ_ctcvcp_4f / THW_ctcvcp_5f, Hincl).

The samples are generated at the inverted top coupling (kHtt = +1 in the model = kappa_t = -1); the
gridpack reweight card stores 69 points in LHEReweightingWeight (see review-notes, Open issues):
    0-15   kHtt +3 ... -3 (without the generated +1), kSM = 1
    16-32  the same scan incl. +1, kSM = 1.5;  33-49  kSM = 0.5;  50-68  CP-mixing scan
In the model's sign convention kHtt = -1 is the SM (kappa_t = +1): index TH_SM_INDEX = 11. Checked on
2018 NanoAODv9 (sigma(SM)/sigma(nominal) = 0.074 tHq, 0.105 tHW).

    th_sm_weight      per-event factor: LHEReweightingWeight[11] for tH datasets, 1 for every other MC
                      dataset; enters the nominal event weight (outvars.nominal_weight_factors)
    normalization     sum of sign(genWeight) * th_sm_weight over ALL generated events (theory_norms:
                      theory_sums[sample][dataset]['all']['nominal'] carries the factor for tH), used by
                      outputToPkl.get_norm instead of sum_signOf_genweights; the cross section stays
                      the SM value, so the tH yield is SM in rate and kinematics
"""
import awkward as ak
import numpy as np

TH_DATASET_PREFIXES = ("THQ_", "THW_")      # dataset names of datasets/th_tzq_definitions.json
TH_SM_INDEX = 11
TH_N_POINTS = 69
# Files dropped from the tH samples (th_sm_weight = 0 for all their events, so they drop out of the yield and of
# the normalization sum): THQ_2017 BC0B1DEE stores 68 points instead of 69 (decided 2026-11-06, review-notes).
# In code rather than only in datasets/THQ.json, since already prepared PBS jobs keep their file lists.
TH_DROPPED_FILES = ("BC0B1DEE-7E4F-0A4F-B385-999056531ADA",)


def is_th(dataset):
    return dataset.startswith(TH_DATASET_PREFIXES)


def th_sm_weight(events, dataset):
    """LHEReweightingWeight[TH_SM_INDEX] for tH datasets, else ones (numpy)."""
    n = len(events)
    if not is_th(dataset):
        return np.ones(n)
    if "LHEReweightingWeight" not in events.fields:
        raise KeyError(f"{dataset}: no LHEReweightingWeight, cannot reweight tH to the SM")
    filename = str(getattr(events, "metadata", {}).get("filename", ""))
    if any(f in filename for f in TH_DROPPED_FILES):
        print(f"--> {dataset}: {filename} dropped (TH_DROPPED_FILES), th_sm_weight = 0")
        return np.zeros(n)
    w = events.LHEReweightingWeight
    counts = ak.to_numpy(ak.num(w))
    if n and not np.all(counts == TH_N_POINTS):
        raise ValueError(f"{dataset}: LHEReweightingWeight has {sorted(set(counts.tolist()))} points, "
                         f"expected {TH_N_POINTS} (reweight card of the ctcvcp samples)")
    return ak.to_numpy(w[:, TH_SM_INDEX]) if n else np.ones(0)
