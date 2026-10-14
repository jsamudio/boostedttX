# Run 3 configuration of the 'trigeff' sample group (samples.GROUPS['trigeff']).
#   pocket-coffea run --cfg inference_conf_trigeff.py -o outputs/nominal/trigeff
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("trigeff")
