import coffea
from coffea.util import load

out = load("/cms/data/store/user/moanwar/ttx_output_hadronic/out_onnx/output_all.coffea")

# Top level keys
print("Top level keys:", out.keys())

# Look at columns structure
print("\nColumns keys:", out['columns'].keys())

# Pick one sample and inspect
sample = list(out['columns'].keys())[0]
print(f"\nFirst sample: {sample}")
print("Datasets:", out['columns'][sample].keys())

dataset = list(out['columns'][sample].keys())[0]
print(f"\nFirst dataset: {dataset}")
print("Categories:", out['columns'][sample][dataset].keys())

cat = list(out['columns'][sample][dataset].keys())[0]
print(f"\nFirst category: {cat}")
print("Collections:", out['columns'][sample][dataset][cat].keys())
