# Run 3 SPANet training inputs (samples.GROUPS['spanet']; samples, cuts and columns in spanet_training.py).
#   pocket-coffea run --cfg inference_conf_spanet.py -o outputs/nominal/spanet
# Then make_spanet_configs.py writes the converter configs (see README.md).
from inference_common import build_config

cfg = build_config("spanet")
