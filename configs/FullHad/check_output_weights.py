import numpy as np
from coffea.util import load

out = load("/cms/data/store/user/moanwar/ttx_output_hadronic/out_onnx/output_all.coffea")

def get_val(col):
    return col.value if hasattr(col, 'value') else np.array(col)

for sample in ['ttHTobb__genMatch', 'QCD_HT']:
    dataset = list(out['columns'][sample].keys())[0]
    cat = 'btag_mask'
    cols = out['columns'][sample][dataset][cat]
    
    print(f"\n=== {sample} / {dataset} ===")
    for key in ['weight', 'events_genWeight', 'events_norm_weight']:
        v = get_val(cols[key])
        print(f"  {key}: mean={np.mean(v):.6f}, min={np.min(v):.6f}, max={np.max(v):.6f}, N={len(v)}")

print("\n=== sum_genweights ===")
for dataset, val in out['sum_genweights'].items():
    print(f"  {dataset}: {val:.2f}")

print("\n=== datasets_metadata ===")
for dataset, meta in out['datasets_metadata'].items():
    print(f"  {dataset}: {meta}")
