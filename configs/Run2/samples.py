"""
Single source of truth for the samples of the analysis.

Pure Python (no coffea / pocket-coffea imports) so that the post-processing scripts
(outputToPkl.py, btag_norm.py, ...) can use it as well as the pocket-coffea configs.

  * DATASET_JSONS    - dataset definition files (datasets/<name>.json)
  * GROUPS           - which samples are processed together by one inference_conf_<group>.py
  * TTBAR_SPLIT /
    GENMATCH_SPLIT   - subsample definitions (turned into Cut objects in inference_common.py)
  * PICKLE_SIGNAL /
    PICKLE_BACKGROUND - which (sub)samples outputToPkl.py exports to the analysis pickles
"""

# ---------------------------------------------------------------- dataset definition files
DATASET_JSONS = [
    "ttHTobb", "ttHToNonbb",
    "TTZToBB", "TTZToQQ", "TTZToLLNuNu",
    "TTbb_Hadronic", "TTbb_2L2Nu", "TTbb_SemiLeptonic",
    "TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonic",
    "QCD", "WJetsToLNu", "DYJetsToLL", "SingleTop", "ttX",
    "DATA_EGamma", "DATA_Muon",
    "ttHSMEFT", "ttZSMEFT", "ttbbSMEFT",
]

# ---------------------------------------------------------------- sample groups
SIGNAL_SAMPLES = ["ttHTobb", "ttHToNonbb", "TTZToQQ", "TTZToLLNuNu", "TTZToBB"]
TTBAR_SAMPLES = ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu",
                 "TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonic"]
MINOR_SAMPLES = ["WJetsToLNu", "DYJetsToLL", "SingleTop", "ttX"]       # add "QCD" here if needed
EFT_SAMPLES = ["ttHSMEFT", "ttZSMEFT", "ttbbSMEFT"]
DATA_SAMPLES = ["DATA_EGamma", "DATA_Muon"]

# One entry per inference_conf_<group>.py. Options:
#   samples   - samples to process
#   is_data   - no MC weights / gen information in the output columns
#   dump_dir  - pocket-coffea 'dump_columns_as_arrays_per_chunk' directory (None = keep in .coffea)
GROUPS = {
    "signal": {"samples": SIGNAL_SAMPLES},
    "bkg":    {"samples": TTBAR_SAMPLES},
    "minor":  {"samples": MINOR_SAMPLES},
    "EFT":    {"samples": EFT_SAMPLES,
               "dump_dir": "/cms/data/store/user/jsamudio/eftchunks/{era}/"},
    "data":   {"samples": DATA_SAMPLES, "is_data": True},
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


def subsample_split(sample):
    """Name of the subsample split used for a sample: 'ttbar', 'genmatch' or None."""
    if sample in TTBAR_SAMPLES:
        return "ttbar"
    if sample in SIGNAL_SAMPLES:
        return "genmatch"
    return None


def subsample_names(sample):
    """Output keys of a sample in the .coffea file, e.g. 'TTToHadronic__tt+LF'."""
    split = subsample_split(sample)
    if split == "ttbar":
        return [f"{sample}__{s}" for s in TTBAR_SPLIT]
    if split == "genmatch":
        return [f"{sample}__{s}" for s in GENMATCH_SPLIT]
    return [sample]


# ---------------------------------------------------------------- sample-specific output columns
def extra_column_group(sample):
    """Which outvars list a sample exports on top of the common columns: 'sig', 'bkg' or None."""
    if sample in SIGNAL_SAMPLES or sample in ("ttHSMEFT", "ttZSMEFT"):
        return "sig"
    if sample in TTBAR_SAMPLES or sample == "ttbbSMEFT":
        return "bkg"
    return None


EFT_COLUMNS = ["EFTfitCoefficients", "EFTfitCoefficientIndex1", "EFTfitCoefficientIndex2", "WCnames"]

# ---------------------------------------------------------------- pickling (outputToPkl.py)
# Both gen-matching categories of every signal sample
PICKLE_SIGNAL = [f"{s}__{m}" for s in ["ttHTobb", "ttHToNonbb", "TTZToBB", "TTZToQQ", "TTZToLLNuNu"]
                 for m in GENMATCH_SPLIT]
# tt+B from the dedicated TTbb samples, tt+LF / tt+C from the inclusive ttbar samples
PICKLE_BACKGROUND = (
    [f"{s}__tt+B" for s in ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu"]]
    + [f"{s}__{c}" for s in ["TTToSemiLeptonic", "TTTo2L2Nu", "TTToHadronic"] for c in ("tt+LF", "tt+C")]
    + ["WJetsToLNu", "DYJetsToLL", "ttX", "SingleTop", "QCD"]
)
