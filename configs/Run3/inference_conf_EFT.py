# Run 3 configuration of the 'EFT' sample group (samples.GROUPS['EFT']).
#   pocket-coffea run --cfg inference_conf_EFT.py -o outputs/nominal/EFT
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("EFT")
