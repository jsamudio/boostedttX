"""
Registry of the SPANet ONNX models (pure Python, shared by the inference configs,
outputToPkl.py and makeDatacard.py).

To switch models, edit ACTIVE_MODELS only. Every active model is evaluated in the
inference step, exported to the pickles (sig/ttbb/ttlf scores -> disc_<model>) and gets
its own datacard.

Model options:
  path            - ONNX file
  jet_features    - AK4 inputs, subset of pt, eta, sin_phi, cos_phi, btagB
  fatjet_features - AK8 inputs, subset of pt, eta, sin_phi, cos_phi, bbtag, xbb_wp_L, xbb_wp_M, xbb_wp_T
  clean_jets      - True if the model was trained on AK4 jets with dR > 0.8 from the
                    leading (highest Xbb) AK8 jet, as done in zhbbworkflow_spanet.py
"""

SPANET_DIR = "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet"

_JET = ["pt", "eta", "sin_phi", "cos_phi", "btagB"]
_FATJET_NOXBB = ["pt", "eta", "sin_phi", "cos_phi"]
_FATJET_XBB = _FATJET_NOXBB + ["bbtag"]
_FATJET_MULTIWP = _FATJET_NOXBB + ["xbb_wp_L", "xbb_wp_M", "xbb_wp_T"]

SPANET_MODELS = {
    "baseline": {
        "path": f"{SPANET_DIR}/run2UL.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_XBB, "clean_jets": False},
    "overlap_NoXbb": {
        "path": f"{SPANET_DIR}/run2UL_NoXbb.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_NOXBB, "clean_jets": False},
    "overlapRemoval_NoXbb": {
        "path": f"{SPANET_DIR}/run2UL_OverlapRemovalNoXbb.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_NOXBB, "clean_jets": True},
    "overlapRemoval_FullXbb": {
        "path": f"{SPANET_DIR}/run2UL_OverlapRemovalFullXbb.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_XBB, "clean_jets": True},
    "overlapRemoval_MultiWP": {
        "path": f"{SPANET_DIR}/run2UL_OverlapRemovalMultiWPXbb.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_MULTIWP, "clean_jets": True},
    "noOverlap_MultiWP": {
        "path": f"{SPANET_DIR}/fullRun2_NoOverlap_MultiWP.onnx",
        "jet_features": _JET, "fatjet_features": _FATJET_MULTIWP, "clean_jets": True},
}

# Models evaluated / exported / turned into datacards
ACTIVE_MODELS = ["noOverlap_MultiWP"]

N_ENCODER_FEATURES = 32


def active_models():
    return {name: SPANET_MODELS[name] for name in ACTIVE_MODELS}


def score_columns(model):
    return [f"sig_score_{model}", f"ttbb_score_{model}", f"ttlf_score_{model}"]


def encoder_columns(model):
    return [f"encoderFeature{i + 1}_{model}" for i in range(N_ENCODER_FEATURES)]
