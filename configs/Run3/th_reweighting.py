"""
SM reweighting of the 2024 tH samples (THQ-4FS-ctcvcp / THW-5FS-ctcvcp), applied by outputToPkl.

The samples are generated at kappa_t = -1 (the model's kHtt = +1) with the Run 2 reweight card (69 points,
SM = kHtt -1, kSM 1 = rwgt_12, index 11). NanoAODv15 does not store these weights (LHEReweightingWeight
is empty in every file), so the SM weight is taken from MiniAOD by extract_th_sm_weights.py (CMSSW) into
one lookup table per dataset, th_sm_weights/<dataset>.npz (key = (luminosityBlock << 32) | event, sorted,
with w and the sign of the event weight; not in git). Without the reweighting the selected tH yield is
1.5-2x too high and up to ~4x in the highest NN bins (Run2/th_reweighting_check.py, 2026-11-27).

Unlike Run 2 the processor does not use it: the inference runs without the tables, and outputToPkl
matches the exported (run, luminosityBlock, event) to the table:

    th_sm_weight      th_sm_weights(): per-event SM weight; in outvars.nominal_weight_factors
    normalization     th_norm_factor() = sum_table(sign) / sum_table(sign * w): with norm = L sigma /
                      sum_signOf_genweights, norm x factor keeps the SM cross section (sum over all
                      generated events, as the Run 2 theory_sums denominator)
"""
import os

import numpy as np

TH_DATASET_PREFIXES = ("THQ_", "THW_")      # dataset names of datasets/TTX.json
TH_SM_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "th_sm_weights")
_TABLES = {}


def is_th(dataset):
    return dataset.startswith(TH_DATASET_PREFIXES)


def event_keys(lumi, event):
    """uint64 key (luminosityBlock << 32) | event; raises if an event number does not fit in 32 bits."""
    lumi, event = np.asarray(lumi, dtype=np.uint64), np.asarray(event, dtype=np.uint64)
    if event.size and event.max() >= np.uint64(1 << 32):
        raise ValueError("event number >= 2^32: the (lumi, event) key needs another layout")
    return (lumi << np.uint64(32)) | event


def load_table(dataset):
    if dataset not in _TABLES:
        path = os.path.join(TH_SM_DIR, f"{dataset}.npz")
        if not os.path.exists(path):
            raise FileNotFoundError(f"{dataset}: no SM weight table {path}. Make it with extract_th_sm_weights.py "
                                    f"(CMSSW, from MiniAOD) and --merge, or pickle with --th-as-is (not for results)")
        t = np.load(path)
        _TABLES[dataset] = {k: t[k] for k in ("key", "w", "sign", "runs")}
    return _TABLES[dataset]


def th_sm_weights(dataset, run, lumi, event):
    """SM weight per event (numpy) from the table; fails if an event is missing."""
    t = load_table(dataset)
    run = np.asarray(run)
    if not np.all(np.isin(run, t["runs"])):
        raise ValueError(f"{dataset}: runs {np.unique(run).tolist()} not all in the table ({t['runs'].tolist()})")
    k = event_keys(lumi, event)
    keys = t["key"]
    idx = np.searchsorted(keys, k)
    found = (idx < len(keys)) & (keys[np.minimum(idx, len(keys) - 1)] == k)
    if not np.all(found):
        raise KeyError(f"{dataset}: {np.count_nonzero(~found)} of {len(k)} events not in the SM weight table "
                       f"(table incomplete? rerun extract_th_sm_weights.py and --merge)")
    return t["w"][idx].astype(np.float64)


def th_norm_factor(dataset):
    """sum(sign) / sum(sign * w) over every generated event of the table."""
    t = load_table(dataset)
    sign = t["sign"].astype(np.float64)
    return float(sign.sum() / np.sum(sign * t["w"].astype(np.float64)))
