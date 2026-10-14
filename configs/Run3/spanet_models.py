"""
Registry of the SPANet ONNX models of the Run 3 analysis (pure Python).

ACTIVE_MODELS is empty until a model trained on the Run 3 object definitions (PUPPI MET,
UParT b tagging, GloParT XbbVsQCD and regressed mass) exists: the inference then runs without
SPANet. Add the new model below and list it in ACTIVE_MODELS.

Model options (inputs as in the Run 2 registry):
  path            - ONNX file
  jet_features    - AK4 inputs, subset of pt, eta, sin_phi, cos_phi, btagB (UParT B)
  fatjet_features - AK8 inputs, subset of pt, eta, sin_phi, cos_phi, bbtag (GloParT XbbVsQCD)
  clean_jets      - True if trained on AK4 jets with dR > 0.8 from the ZH candidate
"""

SPANET_DIR = "/cms/data/jsamudio/boosted/boostedttX/configs"

SPANET_MODELS = {
    # model of the previous 2024 setup (SemiLeptonic_2024): all AK4 jets, GloParT Xbb, PUPPI MET
    "spanet2024": {
        "path": f"{SPANET_DIR}/spanet2024_encoderFeatures.onnx",
        "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
        "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "bbtag"],
        "clean_jets": False},
}

ACTIVE_MODELS = []

N_ENCODER_FEATURES = 32


def active_models():
    return {name: SPANET_MODELS[name] for name in ACTIVE_MODELS}


def score_columns(model):
    return [f"sig_score_{model}", f"ttbb_score_{model}", f"ttlf_score_{model}"]


def encoder_columns(model):
    return [f"encoderFeature{i + 1}_{model}" for i in range(N_ENCODER_FEATURES)]
