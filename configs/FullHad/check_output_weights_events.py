from coffea.util import load
import numpy as np

out = load("/cms/data/store/user/moanwar/ttx_output_hadronic/out_onnx/output_all.coffea")

def get_val(col):
    return col.value if hasattr(col, "value") else np.array(col)

LUMI = 41.5e3  # pb-1

print("=== Quick normalization check ===")
print("Expected = XS * LUMI * (N_selected / N_total)\n")

checks = {
    "ttHTobb__genMatch":  ("ttHTobb_M125_M125_2017",  0.2934,  7825000),
    "QCD_HT":             ("QCD_HT_500to700_500to700_2017", 29990, 36194860),
    "TTZToBB__genMatch":  ("TTZToBB__2017",           0.1157,  7074000),
}

for sample, (dataset, xs, nevents_total) in checks.items():
    cols = out["columns"][sample][dataset]["btag_mask"]
    w = get_val(cols["weight"])
    sum_w = np.sum(w)
    n_sel = len(w)
    
    # What we get from weights
    print(f"--- {sample} ---")
    print(f"  sum(weight)        = {sum_w:.2f}")
    print(f"  N selected         = {n_sel}")
    print(f"  XS={xs} pb, N_total={nevents_total}")
    print(f"  Expected (no genW) = {xs * LUMI * n_sel / nevents_total:.2f}")
    print(f"  Ratio sum_w/expect = {sum_w / (xs * LUMI * n_sel / nevents_total):.4f}")
    print()
