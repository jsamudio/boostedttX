# Run 3 b-tag efficiency maps (samples.GROUPS['btageff']): AK4 jets of the analysis preselection,
# before the b-tag requirement. Then: python3 derive_btag_eff.py
#   pocket-coffea run --cfg inference_conf_btageff.py -o outputs/nominal/btageff
from inference_common import build_config

cfg = build_config("btageff")
