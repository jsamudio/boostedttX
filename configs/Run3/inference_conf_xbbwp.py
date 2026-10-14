# Run 3 Xbb working points (samples.GROUPS['xbbwp']): all fatjets of the analysis preselection, tt+LF
# ttbar (mistag) and ttH(bb) / ttZ(bb) (efficiency). Then: python3 derive_xbb_wps.py
#   pocket-coffea run --cfg inference_conf_xbbwp.py -o outputs/nominal/xbbwp
from inference_common import build_config

cfg = build_config("xbbwp")
