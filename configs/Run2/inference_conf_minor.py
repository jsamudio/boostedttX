# SPANet inference for the 'minor' sample group (samples.GROUPS['minor']).
#   pocket-coffea run --cfg inference_conf_minor.py -o output_minor
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("minor")
run_options = RUN_OPTIONS
