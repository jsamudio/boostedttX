# SPANet inference for the 'EFT' sample group (samples.GROUPS['EFT']).
#   pocket-coffea run --cfg inference_conf_EFT.py -o output_EFT
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("EFT")
run_options = RUN_OPTIONS
