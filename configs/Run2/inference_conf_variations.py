# SPANet inference for the 'variations' sample group (samples.GROUPS['variations']): ttbar / ttbb hdamp and UE samples.
#   pocket-coffea run --cfg inference_conf_variations.py -o output_variations
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("variations")
run_options = RUN_OPTIONS
