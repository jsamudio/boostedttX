# Run 3 configuration of the 'trigeff_mu' sample group (samples.GROUPS['trigeff_mu']): pass 2 of
# TrigEff.ipynb, muon trigger efficiency, analysis muon (electron trigger SF from pass 1 in TRIGGER_SFS).
#   pocket-coffea run --cfg inference_conf_trigeff_mu.py -o outputs/nominal/trigeff_mu
# Eras / samples / jet variation can be set with PC_ERAS, PC_SAMPLES, PC_JEC_VARIATION (see inference_common.py).
from inference_common import build_config

cfg = build_config("trigeff_mu")
