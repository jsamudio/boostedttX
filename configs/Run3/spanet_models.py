"""
Registry of the SPANet ONNX models of the Run 3 analysis (pure Python).

Models trained on the Run 3 object definitions (PUPPI MET, UParT b tagging, GloParT XbbVsQCD and
regressed mass) are listed in ACTIVE_MODELS (empty list: the inference runs without SPANet). Their
ONNX files go to configs/Run3/spanet/models/ (MODEL_DIR) of the checkout that runs the jobs.

Model options (inputs as in the Run 2 registry):
  path            - ONNX file
  jet_features    - AK4 inputs, subset of pt, eta, sin_phi, cos_phi, btagB (UParT B)
  fatjet_features - AK8 inputs, subset of pt, eta, sin_phi, cos_phi, bbtag (GloParT XbbVsQCD),
                    xbb_wp_L, xbb_wp_M, xbb_wp_T (bbtag >= the model's xbb_wps)
  xbb_wps         - {"L": cut, "M": cut, "T": cut} the model was trained with (spanet_training.XBB_WPS);
                    required for xbb_wp_* inputs. Per model, so models with different WPs can run together
  clean_jets      - True if trained on AK4 jets with dR > 0.8 from the ZH candidate
  event_features  - optional event-level inputs (SPANet input "Event"), subset in order of
                    sphericity, aplanarity, centrality, n_jets, n_bjets_<WP> (../common/event_variables.py)
  count_clean_jets - event_features: count only AK4 jets outside the candidate (spanet_training.COUNT_CLEAN_JETS)
The entry for a newly trained model is printed by ../common/make_spanet_configs.py.
"""

import os

SPANET_DIR = "/cms/data/jsamudio/boosted/boostedttX/configs"
# models trained in this setup: configs/Run3/spanet/models/ of the checkout running the jobs
MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spanet", "models")

SPANET_MODELS = {
    # first Run 3 training (2026-11-25), aligned with the Run 2 refactor_noOverlap_MultiWP model;
    # options_run3.json, provisional Xbb WPs (spanet_training.XBB_WPS at training time)
    "run3_refactor": {
        "path": f"{MODEL_DIR}/run3_refactor.onnx",
        "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
        "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "xbb_wp_L", "xbb_wp_M", "xbb_wp_T"],
        "xbb_wps": {"L": 0.25, "M": 0.877, "T": 0.975},
        "event_features": ["sphericity", "aplanarity", "centrality", "n_jets", "n_bjets_L", "n_bjets_M", "n_bjets_T"],
        "count_clean_jets": False,
        "clean_jets": True},
    # model of the previous 2024 setup (SemiLeptonic_2024): all AK4 jets, GloParT Xbb, PUPPI MET
    "spanet2024": {
        "path": f"{SPANET_DIR}/spanet2024_encoderFeatures.onnx",
        "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
        "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "bbtag"],
        "clean_jets": False},
}

ACTIVE_MODELS = ["run3_refactor"]

N_ENCODER_FEATURES = 32

# Classes of the SPANet "EVENT/signal" output (spanet_training.CLASSES):
#   0 ttHTobb, 1 TTZ(->bb), 2 tt+B, 3 tt+C and tt+LF
# sig_score = ttH + ttZ; the single-class scores are stored too
SPANET_CLASSES = {"tth": 0, "ttz": 1, "ttbb": 2, "ttlf": 3}


def active_models():
    return {name: SPANET_MODELS[name] for name in ACTIVE_MODELS}


def score_columns(model):
    return [f"sig_score_{model}", f"ttbb_score_{model}", f"ttlf_score_{model}",
            f"tth_score_{model}", f"ttz_score_{model}"]


def encoder_columns(model):
    return [f"encoderFeature{i + 1}_{model}" for i in range(N_ENCODER_FEATURES)]
