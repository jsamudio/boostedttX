"""
Single source of truth for the samples of the Run 3 analysis (2024 so far).

Pure Python (no coffea / pocket-coffea imports): used by the pocket-coffea configs and by the
post-processing scripts.

  * DATASET_JSONS  - dataset definition files (datasets/<name>.json)
  * GROUPS         - which samples are processed together by one inference_conf_<group>.py
  * PROCESS_LABELS - the 'process' label of each non-ttbar sample (ttbar samples are labelled
                     per event from genTtbarId, see genmatcher.match_gen_tt)
  * TTBAR_SPLIT /
    GENMATCH_SPLIT - subsample definitions (turned into Cut objects in inference_common.py)
"""

import re

import spanet_training

# ---------------------------------------------------------------- dataset definition files
DATASET_JSONS = [
    "ttHTobb", "ttHToNonbb",
    "TTZToQQ", "TTLL", "TTNuNu",
    "TTbb_Hadronic", "TTbb_2L2Nu", "TTbb_SemiLeptonic",
    "TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonic",
    "WJets", "DYJets", "SingleTop", "TTX", "VV", "QCD",
    "DATA_EGamma", "DATA_Muon",
    "ttHSMEFT", "ttbbSMEFT",
]

# ---------------------------------------------------------------- sample groups
# ttZ: there is no dedicated Z->bb sample in Run 3, so Z->bb comes from TTZToQQ (kept as signal);
# Z->ll / Z->nunu from TTLL (two m(ll) bins) and TTNuNu
SIGNAL_SAMPLES = ["ttHTobb", "ttHToNonbb", "TTZToQQ", "TTLL", "TTNuNu"]
TTBAR_SAMPLES = ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu",
                 "TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonic"]

# ttbar / ttbb modelling variations (dedicated samples, same cross section as the nominal sample).
# datasets/<name>.json are built with pocket-coffea build-datasets from
# datasets/ttbar_variations_definitions.json; they are loaded only when processed (inference_common.py).
# No hdamp samples exist for the ttbb in 2024 (DCTR reweighting instead, dctr.py).
TTBAR_VARIATIONS = ["hdampUP", "hdampDOWN", "TuneCP5up", "TuneCP5down"]
TTBB_VARIATIONS = ["TuneCP5up", "TuneCP5down"]
TTBAR_VARIATION_SAMPLES = ([f"{s}_{v}" for s in ["TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu"] for v in TTBAR_VARIATIONS]
                           + [f"{s}_{v}" for s in ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu"] for v in TTBB_VARIATIONS])
MINOR_SAMPLES = ["WJets", "DYJets", "SingleTop", "TTX", "VV"]
QCD_SAMPLES = ["QCD"]
EFT_SAMPLES = ["ttHSMEFT", "ttbbSMEFT"]
DATA_SAMPLES = ["DATA_EGamma", "DATA_Muon"]
# Trigger efficiency (orthogonal method in ttbar e-mu): dileptonic + semileptonic ttbar and both
# data streams (the notebook uses one primary dataset per measurement)
TRIGEFF_SAMPLES = ["TTTo2L2Nu", "TTToSemiLeptonic", "TTbb_2L2Nu", "TTbb_SemiLeptonic"] + DATA_SAMPLES

# Xbb working points (derive_xbb_wps.py): mistag rate on light-flavour fatjets of tt+LF events,
# efficiency on H/Z -> bb fatjets of ttH(bb) / ttZ(bb)
XBBWP_MISTAG_SAMPLES = ["TTToSemiLeptonic", "TTTo2L2Nu"]
XBBWP_SIGNAL_SAMPLES = ["ttHTobb", "TTZToQQ"]
XBBWP_SAMPLES = XBBWP_MISTAG_SAMPLES + XBBWP_SIGNAL_SAMPLES

# b-tag efficiency maps (derive_btag_eff.py), one per group of processes: the jets of each sample
# are corrected with the map of its group. The subsamples are those the analysis takes from the
# sample (tt+B from the 4FS TTbb, tt+LF / tt+C from the 5FS TTTo), so no event enters twice.
BTAG_EFF_GROUPS = {
    "ttLF_C": {"samples": ["TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu"], "subsamples": ["tt+LF", "tt+C"]},
    "ttB":    {"samples": ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu", "ttbbSMEFT"], "subsamples": ["tt+B"]},
    "signal": {"samples": SIGNAL_SAMPLES + ["ttHSMEFT"]},
    "minor":  {"samples": MINOR_SAMPLES + QCD_SAMPLES},
}
# Processed for the maps: the nominal MC of the analysis (QCD uses the 'minor' map)
BTAG_EFF_SAMPLES = SIGNAL_SAMPLES + TTBAR_SAMPLES + MINOR_SAMPLES


def btag_eff_group(sample, subsample=None):
    """Efficiency-map group of a sample (variation samples: that of their nominal sample). With a
    subsample: None if the analysis does not use that subsample of the sample."""
    nominal = re.sub(r"_(hdampUP|hdampDOWN|TuneCP5up|TuneCP5down)$", "", sample)
    for group, opts in BTAG_EFF_GROUPS.items():
        if nominal in opts["samples"]:
            if subsample is not None and subsample not in opts.get("subsamples", [subsample]):
                return None
            return group
    raise ValueError(f"Sample '{sample}' is in no b-tag efficiency group (samples.BTAG_EFF_GROUPS)")


# One entry per inference_conf_<group>.py. Options:
#   samples   - samples to process
#   is_data   - no MC weights / gen information in the output columns
#   selection - "analysis" (default), "trigeff" (e-mu ttbar, no trigger requirement) or "spanet"
#               (SPANet training inputs, samples and options in spanet_training.py) or "xbbwp"
#               (all fatjets of the preselection, for the Xbb working points) or "btageff"
#               (AK4 jets of the preselection, b-tag efficiency histograms)
#   dump_dir  - pocket-coffea 'dump_columns_as_arrays_per_chunk' directory (None = keep in .coffea)
GROUPS = {
    "signal":  {"samples": SIGNAL_SAMPLES},
    "bkg":     {"samples": TTBAR_SAMPLES},
    "minor":   {"samples": MINOR_SAMPLES},
    "variations": {"samples": TTBAR_VARIATION_SAMPLES},
    "qcd":     {"samples": QCD_SAMPLES},
    "EFT":     {"samples": EFT_SAMPLES},
    "data":    {"samples": DATA_SAMPLES, "is_data": True},
    "trigeff": {"samples": TRIGEFF_SAMPLES, "selection": "trigeff"},
    "spanet":  {"samples": spanet_training.sample_names(), "selection": "spanet"},
    "xbbwp":   {"samples": XBBWP_SAMPLES, "selection": "xbbwp"},
    "btageff": {"samples": BTAG_EFF_SAMPLES, "selection": "btageff"},
}

# ---------------------------------------------------------------- process labels
PROCESS_LABELS = {
    "ttHTobb": "ttH", "ttHToNonbb": "ttH", "ttHSMEFT": "ttH",
    "TTZToQQ": "ttZ", "TTLL": "ttZ", "TTNuNu": "ttZ",
    "WJets": "VJets", "DYJets": "VJets",
    "SingleTop": "SingleTop", "TTX": "ttX", "VV": "VV", "QCD": "QCD",
    "DATA_EGamma": "data_obs", "DATA_Muon": "data_obs",
}

# ---------------------------------------------------------------- subsamples
# genTtbarId % 100 categories (https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher)
TTBAR_SPLIT = {
    "tt+LF": [0],
    "tt+C":  [41, 42, 43, 44, 45, 46],
    "tt+B":  [51, 52, 53, 54, 55, 56],
}
# Signal split on matchedGen_ZHbb_bb (see genmatcher.match_gen_sig)
GENMATCH_SPLIT = ["genMatch", "non_genMatch"]


def subsample_names(sample):
    """Output keys of a sample in the .coffea file, e.g. 'TTToHadronic__tt+LF'."""
    if sample in TTBAR_SAMPLES or sample in TTBAR_VARIATION_SAMPLES:
        return [f"{sample}__{s}" for s in TTBAR_SPLIT]
    if sample in SIGNAL_SAMPLES:
        return [f"{sample}__{s}" for s in GENMATCH_SPLIT]
    return [sample]


# ---------------------------------------------------------------- sample-specific output columns
def extra_column_group(sample):
    """Which outvars list a sample exports on top of the common columns: 'sig', 'bkg' or None."""
    if sample in SIGNAL_SAMPLES or sample == "ttHSMEFT":
        return "sig"
    if sample in TTBAR_SAMPLES or sample in TTBAR_VARIATION_SAMPLES or sample == "ttbbSMEFT":
        return "bkg"
    return None


EFT_COLUMNS = ["EFTfitCoefficients", "EFTfitCoefficientIndex1", "EFTfitCoefficientIndex2", "WCnames"]
