#!/usr/bin/env python3
"""
merge_coffea.py — merge per-sample coffea outputs into one file.
Called by run_inference.sh after all samples for a year/type complete.

Usage:
    python3 merge_coffea.py --outdir /path/to/outdir --samples ttHTobb TTZToBB QCD_HT ...
"""

import os
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir",  required=True,
                        help="Directory containing _<sample>/ subdirs")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="Sample names to merge")
    args = parser.parse_args()

    try:
        from coffea.util import load, save
        from coffea.processor import accumulate
    except ImportError:
        print("[ERROR] coffea not available")
        sys.exit(1)

    final_output = os.path.join(args.outdir, "output_all.coffea")
    outputs  = []
    missing  = []

    for sample in args.samples:
        coffea_file = os.path.join(args.outdir, f"_{sample}", "output_all.coffea")
        if os.path.exists(coffea_file):
            print(f"  Loading: {coffea_file}")
            try:
                out = load(coffea_file)
                outputs.append(out)
            except Exception as e:
                print(f"  [ERROR] Failed to load {coffea_file}: {e}")
                missing.append(sample)
        else:
            print(f"  [MISSING] {coffea_file}")
            missing.append(sample)

    if not outputs:
        print("[ERROR] No outputs to merge!")
        sys.exit(1)

    print(f"\nMerging {len(outputs)} coffea files...")
    merged = accumulate(outputs)
    save(merged, final_output)
    print(f"[OK] Saved: {final_output}")

    # Print dataset summary
    if "sum_signOf_genweights" in merged:
        datasets = sorted(merged["sum_signOf_genweights"].keys())
        print(f"\nDatasets in merged output ({len(datasets)}):")
        for d in datasets:
            print(f"  {d}")

    if missing:
        print(f"\n[WARNING] Missing samples not included in merge: {missing}")
    else:
        print(f"\n[SUCCESS] All {len(outputs)} samples merged successfully")


if __name__ == "__main__":
    main()
