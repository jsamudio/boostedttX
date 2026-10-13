# SPANet inference for the 'signal' sample group (samples.GROUPS['signal']).
#   pocket-coffea run --cfg inference_conf_signal.py -o output_signal
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config, RUN_OPTIONS

cfg = build_config("signal")
run_options = RUN_OPTIONS
