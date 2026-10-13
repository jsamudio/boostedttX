# SPANet inference for the 'bkg' sample group (samples.GROUPS['bkg']).
#   pocket-coffea run --cfg inference_conf_bkg.py -o output_bkg
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("bkg")
run_options = RUN_OPTIONS
