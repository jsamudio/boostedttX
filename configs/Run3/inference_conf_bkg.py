# Run 3 configuration of the 'bkg' sample group (samples.GROUPS['bkg']).
#   pocket-coffea run --cfg inference_conf_bkg.py -o outputs/nominal/bkg
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("bkg")
