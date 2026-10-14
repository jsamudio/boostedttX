# Run 3 configuration of the 'data' sample group (samples.GROUPS['data']).
#   pocket-coffea run --cfg inference_conf_data.py -o outputs/nominal/data
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("data")
