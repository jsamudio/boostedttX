# Run 3 configuration of the 'trigeff_ele' sample group (samples.GROUPS['trigeff_ele']): pass 1 of
# TrigEff.ipynb, electron trigger efficiency, reference muon = MUO trigger-SF denominator.
#   pocket-coffea run --cfg inference_conf_trigeff_ele.py -o outputs/nominal/trigeff_ele
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("trigeff_ele")
