# Run 3 configuration of the 'qcd' sample group (samples.GROUPS['qcd']).
#   pocket-coffea run --cfg inference_conf_qcd.py -o outputs/nominal/qcd
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("qcd")
