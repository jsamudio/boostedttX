# Run 3 configuration of the 'signal' sample group (samples.GROUPS['signal']).
#   pocket-coffea run --cfg inference_conf_signal.py -o outputs/nominal/signal
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("signal")
