# Run 3 configuration of the 'variations' sample group (samples.GROUPS['variations']): ttbar hdamp / UE and ttbb UE samples.
#   pocket-coffea run --cfg inference_conf_variations.py -o outputs/nominal/variations
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("variations")
