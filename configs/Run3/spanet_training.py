"""
SPANet training inputs of the Run 3 analysis: the ONE place to edit when the training inputs change.

Pure Python (no coffea imports). Read by
  * inference_common.build_config("spanet") / inference_conf_spanet.py  (pocket-coffea columns, cuts,
    subsamples, weights, gen matching options)
  * make_spanet_configs.py  (YAML configs of the converters in ../spanet_inputs: coffea -> parquet -> h5)

Workflow (see README.md, "SPANet training inputs"):
    python3 run_pipeline.py inference --backend pbs --groups spanet        # 1. ntuples (.coffea)
    python3 make_spanet_configs.py                                         # 2. converter YAMLs
    cd ../spanet_inputs
    python3 coffea_to_parquet_script.py -c ../Run3/spanet/features_parquet.yaml \\
        -i ../Run3/outputs/nominal/spanet/output_all.coffea -o <dir>/run3.parquet --cat spanet
    python3 parquet_to_h5.py -c ../Run3/spanet/features_h5.yaml -i <dir>/run3_*.parquet -o <dir>/run3.h5

Gen matching (MC, spanet_training_workflow.py): hard-process partons from GenPart (ancestry, not the
LHE ordering, so it works for Powheg ttH and MadGraph ttZ+jet alike) with the Run 2 provenance codes
  1 = H/Z -> qq (and the extra g -> bb in the 4F ttbb sample, Run 2 convention)
  2 = b from the hadronic top    3 = b from the leptonic top
  4 = additional radiation       5 = quarks from the hadronic W      -1 = unmatched
AK4 jets get `matched` / `prov` / `dRMatchedJet` from a dR match to these partons; the fatjets get
`matched` / `prov` from a dR match to the hard-process H/Z (-> qq).
"""

# ---------------------------------------------------------------- samples and classes
# Training samples (must be in samples.DATASET_JSONS). Options:
#   subsamples         - None (whole sample), "ttbar" (genTtbarId split: tt+LF / tt+C / tt+B) or
#                        {name: event flag} with boolean event columns set by the workflow
#                        (genmatcher.match_gen_sig: Zbb, Hbb, Zqq, ...)
#   extra_b_as_boson   - label hard-process b quarks not from t/W/Z/H with provenance 1 (4F ttbb)
SAMPLES = {
    "ttHTobb":           {"subsamples": None},
    "TTZToQQ":           {"subsamples": {"Zbb": "Zbb"}},        # no dedicated TTZ->bb sample in Run 3
    "TTbb_SemiLeptonic": {"subsamples": "ttbar", "extra_b_as_boson": True},
    "TTToSemiLeptonic":  {"subsamples": "ttbar"},
}

# Classification targets: <sample>__<subsample> (or sample) substring -> (class, one-hot index).
# The converter writes one parquet per (sub)sample, named <output>_<sample>__<subsample>.parquet; each
# key must match exactly one of those names.
CLASSES = {
    "ttHTobb":      ("tthbb", 1),
    "TTZToQQ__Zbb": ("ttzbb", 2),
    "tt+B":         ("ttbb", 3),
    "tt+C":         ("ttlf", 4),
    "tt+LF":        ("ttlf", 4),
}

# ---------------------------------------------------------------- event selection
# The analysis preselection (trigger, 1 lepton, >= 5 AK4, >= 1 AK8, PUPPI MET, low-mass vetoes) and,
# on top, the category below. None disables a requirement.
CATEGORY = "spanet"
SELECTION = {
    "min_b_jets":  2,           # nbJetGood >= (UParT medium), the analysis btag_mask
    "min_ak4_jets": 5,          # n_ak4jets >= (after the optional cleaning below)
    "n_b_outZH":   2,           # b jets outside the candidate cone ==
    "zh_mass":     (50, 200),   # GloParT regressed mass window of the candidate [GeV]
    "min_xbb":     None,        # candidate GloParT XbbVsQCD >= (no Run 3 WP yet)
}

# AK4 jets used by SPANet: False = all JetGood, True = only jets with dR > 0.8 from the ZH candidate.
# Must agree with the model's "clean_jets" in spanet_models.py.
CLEAN_JETS = False

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
# "events" holds per-event columns. Fields added by the workflow: JetGood.btag_{L,M,T,XT,XXT},
# JetGood/FatJetSorted.{matched,prov,dRMatchedJet}, LeptonGood.is_electron, events.n_*.
EXPORT = {
    "JetGood": {"jagged": True, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass",
        "btagB": "btagB", "btag_L": "btag_L", "btag_M": "btag_M", "btag_T": "btag_T",
        "hadronFlavour": "hadronFlavour",
        "matched": "matched", "prov": "prov", "dRMatchedJet": "dRMatchedJet"}},
    "FatJetSorted": {"max": 1, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "mass": "mass", "bbtag": "btagBB",
        "matched": "matched", "prov": "prov", "dRMatchedJet": "dRMatchedJet"}},
    "LeptonGood": {"max": 1, "fields": {
        "pt": "pt", "eta": "eta", "phi": "phi", "pdgId": "pdgId", "is_electron": "is_electron"}},
    "PuppiMET": {"fields": {"pt": "pt", "phi": "phi"}},
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
        "nBosonPartonsInFatJet": "nBosonPartonsInFatJet"}},
}

# Matched collections for the parquet converter: entries with pt == -999 become None
MATCHED_COLLECTIONS = {"Parton": "PartonMatched"}

# ---------------------------------------------------------------- SPANet inputs (h5)
# {SPANet input: (exported collection, [features])}. Besides the exported feature names the h5
# converter provides MASK (pt != 0), sin_phi, cos_phi and is_electron (from pdgId).
# The inference (spanet_inference.onnx_inference) feeds: Jet [pt, eta, sin_phi, cos_phi, btagB] +
# mask, Lepton [pt, eta, sin_phi, cos_phi], Met [pt, sin_phi, cos_phi], FatJet [pt, eta, sin_phi,
# cos_phi, bbtag]; a model trained on other inputs needs the inference extended accordingly
# (make_spanet_configs.py warns).
SPANET_INPUTS = {
    "Jet":    ("JetGood",      ["MASK", "pt", "eta", "sin_phi", "cos_phi", "btagB"]),
    "Lepton": ("LeptonGood",   ["pt", "eta", "sin_phi", "cos_phi"]),
    "Met":    ("PuppiMET",     ["pt", "sin_phi", "cos_phi"]),
    "FatJet": ("FatJetSorted", ["pt", "eta", "sin_phi", "cos_phi", "bbtag"]),
}
PARTICLES = ["t1", "t2", "h"]       # assignment targets (h5.py: t1 = hadronic top, t2 = leptonic b, h = prov 1)
FRAC_TRAIN = 0.8


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
