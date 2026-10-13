import awkward as ak
import numpy as np
import os
import json

# ── Lumi map (pb-1) — must be at module level ─────────────────────────────────
LUMI_MAP = {
    "2016_PreVFP":  19500.,
    "2016_PostVFP": 16800.,
    "2017":         41500.,
    "2018":         59740.,
}


# ── TTbb sum(|genWeight|) lookup ──────────────────────────────────────────────
_TTBB_SUMW_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "ttbb_sumw.json")
_TTBB_SUMW = {}

if os.path.exists(_TTBB_SUMW_FILE):
    with open(_TTBB_SUMW_FILE) as f:
        _data = json.load(f)
    for key, v in _data.items():
        _TTBB_SUMW[v["dataset_key"]] = v["sum_abs_genw"]
    print(f"[weight_handler] Loaded TTbb sum_abs_genw for "
          f"{len(_TTBB_SUMW)} datasets")
else:
    print(f"[weight_handler] WARNING: ttbb_sumw.json not found — "
          f"TTbb weights will be wrong!")
    print(f"  Run: python3 compute_ttbb_sumw.py --all")


def calc_weight(events, output, dataset, params):
    if events.metadata['isMC'] == 'True':
        year    = events.metadata['year']
        lumi    = LUMI_MAP.get(year, params.sample_params['lumi']['lumi'] * 1000)
        xs_lumi = float(events.metadata['xsec']) * lumi

        if dataset in _TTBB_SUMW:
            # TTbb: per-event norm_weight so PocketCoffea formula gives correct result
            # w = norm_weight / sum_signOf_genweights × sign(genW)
            #   = xs_lumi × |genW| / sum_abs / sum_signOf_genweights × sign(genW)
            #   = xs_lumi × genW / sum_abs / sum_signOf_genweights × sum_signOf_genweights
            # Wait — this still has sum_signOf_genweights in denominator from plotting
            # See note below
            sum_abs = _TTBB_SUMW[dataset]
            genw    = ak.to_numpy(events.genWeight).astype(np.float64)
            norm_w  = xs_lumi * np.abs(genw) / sum_abs
            events['norm_weight'] = ak.from_numpy(norm_w.astype(np.float32))
        else:
            # Standard samples: constant norm_weight = XS × lumi
            events['norm_weight'] = xs_lumi
    else:
        events['norm_weight'] = 1
