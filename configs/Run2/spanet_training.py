"""
SPANet training inputs of the Run 2 analysis: the ONE place to edit when the training inputs change.

Pure Python (no coffea imports). Read by
  * inference_common.build_config("spanet") / inference_conf_spanet.py  (pocket-coffea columns, cuts,
    subsamples, weights, gen matching options)
  * ../common/make_spanet_configs.py  (YAML configs of the converters in ../spanet_inputs: coffea ->
    parquet -> h5, and the INPUTS section of the SPANet event-info file)

Workflow (see README.md, "SPANet training inputs"):
    python3 run_pipeline.py inference --backend pbs --groups spanet --eras <eras>   # 1. ntuples (.coffea)
    python3 ../common/make_spanet_configs.py                                        # 2. converter YAMLs
    cd ../spanet_inputs
    python3 coffea_to_parquet_script.py -c ../Run2/spanet/features_parquet.yaml \\
        -i ../Run2/outputs/nominal/spanet/output_all.coffea -o <dir>/run2.parquet --cat spanet
    python3 parquet_to_h5.py -c ../Run2/spanet/features_h5.yaml -i <dir>/run2_*.parquet -o <dir>/run2.h5

Same objects and preselection as the analysis (spanet_inference.ZHbbSpanetProcessor). Gen matching
(../common/parton_matching.py, shared with Run 3): hard-process partons from the GenPart ancestry with
the provenance codes
  1 = H/Z -> qq (and the extra g -> bb in the 4F ttbb sample)
  2 = b from the hadronic top    3 = b from the leptonic top
  4 = additional radiation       5 = quarks from the hadronic W      -1 = unmatched
"""

# ---------------------------------------------------------------- samples and classes
# Training samples (must be in samples.DATASET_JSONS). Options:
#   subsamples         - None (whole sample), "ttbar" (genTtbarId split: tt+LF / tt+C / tt+B) or
#                        {name: event flag} with boolean event columns set by the workflow
#                        (genmatcher.match_gen_sig: Zbb, Hbb, Zqq, ...)
#   extra_b_as_boson   - label hard-process b quarks not from t/W/Z/H with provenance 1 (4F ttbb)
SAMPLES = {
    "ttHTobb":           {"subsamples": {"genMatch": "matchedGen_ZHbb_bb"}},
    "TTZToBB":           {"subsamples": {"genMatch": "matchedGen_ZHbb_bb"}},
    "TTbb_SemiLeptonic": {"subsamples": "ttbar", "extra_b_as_boson": True},
    "TTToSemiLeptonic":  {"subsamples": "ttbar"},
}

# Classification targets: <sample>__<subsample> (or sample) substring -> (class, index). The index is
# the SPANet output position (0..N-1), which spanet_models.SPANET_CLASSES must follow.
CLASSES = {
    "ttHTobb": ("tthbb", 0),
    "TTZToBB": ("ttzbb", 1),
    "tt+B":    ("ttbb", 2),
    "tt+C":    ("ttlf", 3),
    "tt+LF":   ("ttlf", 3),
}

# ---------------------------------------------------------------- event selection
# The analysis preselection (triggers, 1 lepton, >= 5 AK4, >= 1 AK8, MET, low-mass vetoes, HEM) and,
# on top, the category below (the old training's btag_mask + network_cut). None disables a requirement.
CATEGORY = "spanet"
SELECTION = {
    "min_b_jets":  2,           # nbJetGood >= (DeepJet medium), the analysis btag_mask
    "min_ak4_jets": 5,          # n_ak4jets >= (all AK4 jets, as in the analysis)
    "n_b_outZH":   2,           # b jets outside the candidate cone ==
    "zh_mass":     (50, 200),   # ParticleNet regressed mass window of the candidate [GeV]
    "xbb_wp":      "L",         # candidate XbbVsQCD >= this working point of the era (analysis_years)
}

# AK4 jets used by SPANet: False = all JetGood, True = only jets with dR > 0.8 from the ZH candidate.
# Must agree with the model's "clean_jets" in spanet_models.py.
CLEAN_JETS = True

# Event-level inputs (../common/event_variables.py): sphericity, aplanarity, centrality over the
# lepton, the candidate fat jet and the AK4 jets outside it; jet / b-jet counts.
# COUNT_CLEAN_JETS: count only AK4 jets outside the candidate (False = all AK4 jets).
# Must agree with the model's "event_features" / "count_clean_jets" in spanet_models.py.
BTAG_COUNT_WPS = ["L", "M", "H"]        # DeepJet working points (pocket-coffea: L, M, H)
COUNT_CLEAN_JETS = False
EVENT_FEATURES = ["sphericity", "aplanarity", "centrality", "n_jets"] + [f"n_bjets_{wp}" for wp in BTAG_COUNT_WPS]

# ---------------------------------------------------------------- gen matching
MATCHING = {
    "jet_dr":     0.3,              # parton - AK4 jet
    "fatjet_dr":  0.6,              # H/Z - AK8 jet
    "parton_copy": "isFirstCopy",   # GenPart copy used for the partons ("isFirstCopy" ~ LHE level,
                                    #  "isLastCopyBeforeFSR" = after ISR, before FSR)
}

# Event weight of the "weight" column (pocket-coffea weights); the converter divides it by the sum
# of genweights, so it is normalized to xsec x lumi per sample
WEIGHTS = ["genWeight", "lumi", "XS", "pileup"]

# ---------------------------------------------------------------- exported collections
# {collection: {"fields": {feature name: field}, "jagged": bool, "max": n}}
#   feature name - name in the parquet / h5 files (and the SPANet input names)
#   field        - field of the collection in the workflow
#   jagged       - keep as a variable-length list (awkward collection in the parquet)
#   max          - export only the first n objects (1 = one value per event, as a flat column)
# "events" holds per-event columns. Fields added by the workflow: JetGood.btag_{L,M,H},
# JetGood/FatJetSorted.{matched,prov,dRMatchedJet}, LeptonGood.is_electron, FatJetSorted.xbb_wp_*,
# events.n_*, the EVENT_FEATURES.
EXPORT = {
    "JetGood": {"jagged": True, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass",
        "btagB": "btagB", "btag_L": "btag_L", "btag_M": "btag_M", "btag_H": "btag_H",
        "hadronFlavour": "hadronFlavour",
        "matched": "matched", "prov": "prov", "dRMatchedJet": "dRMatchedJet"}},
    "FatJetSorted": {"max": 1, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass", "bbtag": "xbbVsQCD",
        "xbb_wp_L": "xbb_wp_L", "xbb_wp_M": "xbb_wp_M", "xbb_wp_T": "xbb_wp_T",
        "matched": "matched", "prov": "prov", "dRMatchedJet": "dRMatchedJet"}},
    "LeptonGood": {"max": 1, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "pdgId": "pdgId", "is_electron": "is_electron"}},
    "MET": {"fields": {"pt": "pt", "phi": "phi"}},
    "Parton": {"jagged": True, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass", "pdgId": "pdgId", "prov": "provenance"}},
    # partons in the AK4-jet ordering (None -> -999 where the jet is unmatched)
    "PartonMatched": {"jagged": True, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass", "pdgId": "pdgId", "prov": "provenance"}},
    "LeptonParton": {"jagged": True, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass", "pdgId": "pdgId"}},
    "Generator": {"fields": {"x1": "x1", "x2": "x2", "pdgid1": "id1", "pdgid2": "id2"}},
    "events": {"fields": {
        "n_ak4jets": "n_ak4jets", "n_b_inZH": "n_b_inZH", "n_b_outZH": "n_b_outZH", "ZH_M": "ZH_M",
        "nParton": "nParton", "nPartonMatched": "nPartonMatched",
        "nBosonPartonsInFatJet": "nBosonPartonsInFatJet",
        **{f: f for f in EVENT_FEATURES}}},
}

# Matched collections for the parquet converter: entries with pt == -999 become None
MATCHED_COLLECTIONS = {"Parton": "PartonMatched"}

# ---------------------------------------------------------------- SPANet inputs (h5)
# {SPANet input: (exported collection, [features])}. Besides the exported feature names the h5
# converter provides MASK (pt != 0), sin_phi, cos_phi and is_electron (from pdgId).
SPANET_INPUTS = {
    "Jet":    ("JetGood",      ["MASK", "pt", "eta", "sin_phi", "cos_phi", "btagB"]),
    "Lepton": ("LeptonGood",   ["pt", "eta", "sin_phi", "cos_phi"]),
    "Met":    ("MET",          ["pt", "sin_phi", "cos_phi"]),
    "FatJet": ("FatJetSorted", ["pt", "eta", "sin_phi", "cos_phi", "xbb_wp_L", "xbb_wp_M", "xbb_wp_T"]),
    "Event":  ("events",       list(EVENT_FEATURES)),
}
PARTICLES = ["t1", "t2", "h"]       # assignment targets (h5.py: t1 = hadronic top, t2 = leptonic b, h = prov 1)
FRAC_TRAIN = 0.8

# What the inference (spanet_inference.onnx_inference) can feed: Jet / FatJet / Event any subset (the
# model's feature lists in spanet_models.py), Lepton and Met fixed
INFERENCE_FEATURES = {
    "Jet": ["MASK", "pt", "eta", "sin_phi", "cos_phi", "btagB"],
    "FatJet": ["pt", "eta", "sin_phi", "cos_phi", "bbtag", "xbb_wp_L", "xbb_wp_M", "xbb_wp_T"],
    "Event": list(EVENT_FEATURES),
    "Lepton": ["pt", "eta", "sin_phi", "cos_phi"],
    "Met": ["pt", "sin_phi", "cos_phi"],
}
INFERENCE_FIXED = {"Lepton", "Met"}


# ---------------------------------------------------------------- helpers
def sample_names():
    return list(SAMPLES)


def exported_fields(collection):
    """Fields of the workflow collection to export (deduplicated, in order)."""
    return list(dict.fromkeys(EXPORT[collection]["fields"].values()))


def validate():
    """Basic consistency checks (raises ValueError)."""
    for name, (collection, features) in SPANET_INPUTS.items():
        if collection not in EXPORT:
            raise ValueError(f"SPANET_INPUTS[{name}] uses {collection}, which is not in EXPORT")
        derived = {"MASK": "pt", "sin_phi": "phi", "cos_phi": "phi", "is_electron": "pdgId"}
        for feature in features:
            needed = derived.get(feature, feature)
            if needed not in EXPORT[collection]["fields"]:
                raise ValueError(f"SPANet input {name}.{feature} needs feature '{needed}' of {collection} in EXPORT")
    for coll, matched in MATCHED_COLLECTIONS.items():
        for c in (coll, matched):
            if c not in EXPORT or not EXPORT[c].get("jagged"):
                raise ValueError(f"MATCHED_COLLECTIONS: {c} must be a jagged collection in EXPORT")
    for sample, opts in SAMPLES.items():
        subs = opts.get("subsamples")
        if not (subs is None or subs == "ttbar" or isinstance(subs, dict)):
            raise ValueError(f"SAMPLES[{sample}]['subsamples'] must be None, 'ttbar' or a dict")
    for collection, spec in EXPORT.items():
        if spec.get("jagged") and spec.get("max") == 1:
            raise ValueError(f"EXPORT[{collection}]: max = 1 is a flat column, drop 'jagged'")
