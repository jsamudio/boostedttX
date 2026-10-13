# SPANet inference for the 'data' sample group (samples.GROUPS['data']).
#   pocket-coffea run --cfg inference_conf_data.py -o output_data
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("data")
run_options = RUN_OPTIONS
