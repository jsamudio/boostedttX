# Run 3 configuration of the 'minor' sample group (samples.GROUPS['minor']).
#   pocket-coffea run --cfg inference_conf_minor.py -o outputs/nominal/minor
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("minor")
