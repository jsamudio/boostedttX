#!/usr/bin/env python3
"""
compute_ttbb_sumw.py
Computes sum(|genWeight|) over ALL events for TTbb samples from NanoAOD files.
This is needed for correct normalisation of TTbb 4FS samples with negative weights.

Usage:
    python3 compute_ttbb_sumw.py --year 2017
    python3 compute_ttbb_sumw.py --all

Output:
    ttbb_sumw.json — dict of {dataset_key: sum_abs_genweights}
"""

import json
import argparse
import numpy as np
import os

LOCALDIR = os.path.dirname(os.path.abspath(__file__))

TTBB_SAMPLES = {
    "TTbb_Hadronic":     "datasets/TTbb_Hadronic.json",
    "TTbb_SemiLeptonic": "datasets/TTbb_SemiLeptonic.json",
    "TTbb_2L2Nu":        "datasets/TTbb_2L2Nu.json",
}

LUMI_MAP = {
    "2016APV": 19500.,
    "2016":    16800.,
    "2017":    41500.,
    "2018":    59740.,
}


def get_file_path(file_entry):
    """
    Extract file path from JSON entry.
    Handles both plain string and dict formats:
      "file:/path/to/file.root"  → "/path/to/file.root"
      {"path": "file:/path/..."}  → "/path/to/file.root"
    """
    if isinstance(file_entry, dict):
        path = file_entry.get("path", file_entry.get("file", ""))
    else:
        path = str(file_entry)
    # Strip file:// or file: prefix
    path = path.replace("file://", "").replace("file:", "")
    return path


def compute_sum_abs_genw(files, max_files=None):
    """Compute sum(|genWeight|) over all events in a list of NanoAOD files."""
    try:
        import uproot
    except ImportError:
        print("  [ERROR] uproot not available. Install with: pip install uproot")
        return None

    total    = 0.0
    n_events = 0
    files_to_process = files[:max_files] if max_files else files

    for i, f in enumerate(files_to_process):
        fpath = get_file_path(f)
        if not os.path.exists(fpath):
            print(f"  [SKIP] File not found: {fpath}")
            continue
        try:
            with uproot.open(f"{fpath}:Events") as tree:
                genw      = tree["genWeight"].array(library="np")
                total    += np.sum(np.abs(genw))
                n_events += len(genw)
                if (i + 1) % 10 == 0 or i == 0:
                    print(f"    [{i+1}/{len(files_to_process)}] "
                          f"n_events={n_events:,}  sum_abs={total:.4e}")
        except Exception as e:
            print(f"  [ERROR] {fpath}: {e}")
            continue

    return total, n_events


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year",      nargs="+",
                        choices=["2016APV", "2016", "2017", "2018"],
                        default=["2017"])
    parser.add_argument("--all",       action="store_true",
                        help="Process all years")
    parser.add_argument("--output",    default="ttbb_sumw.json")
    parser.add_argument("--max-files", type=int, default=None,
                        help="Max files per dataset (for testing, e.g. --max-files 5)")
    args = parser.parse_args()

    years = ["2016APV", "2016", "2017", "2018"] if args.all else args.year

    # Load existing results
    if os.path.exists(args.output):
        with open(args.output) as f:
            results = json.load(f)
        print(f"Loaded existing results from {args.output} "
              f"({len(results)} entries)")
    else:
        results = {}

    for year in years:
        print(f"\n{'='*60}")
        print(f"Year: {year}")
        print(f"{'='*60}")

        for sample_name, json_rel in TTBB_SAMPLES.items():
            json_path = os.path.join(LOCALDIR, json_rel)
            if not os.path.exists(json_path):
                print(f"  [SKIP] JSON not found: {json_path}")
                continue

            d = json.load(open(json_path))

            # Find dataset key for this year
            dataset_key = None
            for k in d:
                if year in k:
                    dataset_key = k
                    break
            if dataset_key is None:
                print(f"  [SKIP] No {year} entry in {json_rel}")
                continue

            result_key = f"{sample_name}_{year}"
            if result_key in results:
                existing = results[result_key]["sum_abs_genw"]
                print(f"  [SKIP] {result_key} already computed: "
                      f"sum_abs={existing:.4e}")
                continue

            xsec      = float(d[dataset_key]["metadata"]["xsec"])
            nevents   = int(d[dataset_key]["metadata"]["nevents"])
            lumi      = LUMI_MAP[year]
            files     = d[dataset_key]["files"]

            print(f"\n  {sample_name} ({year}):")
            print(f"    dataset_key = {dataset_key}")
            print(f"    xsec        = {xsec} pb")
            print(f"    nevents     = {nevents:,}")
            print(f"    lumi        = {lumi} pb-1")
            print(f"    XS*lumi     = {xsec*lumi:.1f}")
            print(f"    n_files     = {len(files)}")

            # Show first file path for debugging
            first_path = get_file_path(files[0]) if files else "N/A"
            print(f"    first_file  = {first_path}")
            if not os.path.exists(first_path):
                print(f"    [WARNING] First file not found on disk!")

            result = compute_sum_abs_genw(files, max_files=args.max_files)
            if result is None:
                continue

            sum_abs, n_proc = result
            if n_proc == 0:
                print(f"  [ERROR] No events processed for {result_key}")
                continue

            mean_abs = sum_abs / n_proc

            print(f"\n    === Results for {result_key} ===")
            print(f"    n_events_processed = {n_proc:,}")
            print(f"    sum_abs_genw       = {sum_abs:.6e}")
            print(f"    mean_abs_genw      = {mean_abs:.4f}")
            print(f"    norm_weight (mean) = {xsec*lumi/sum_abs:.6e}")
            print(f"    expected yield     ≈ {xsec*lumi:.1f} events @ {lumi/1000:.1f} fb-1")

            results[result_key] = {
                "sample":        sample_name,
                "year":          year,
                "dataset_key":   dataset_key,
                "xsec":          xsec,
                "lumi":          lumi,
                "xs_lumi":       xsec * lumi,
                "nevents_json":  nevents,
                "n_processed":   n_proc,
                "sum_abs_genw":  sum_abs,
                "mean_abs_genw": mean_abs,
            }

            with open(args.output, "w") as f:
                json.dump(results, f, indent=2)
            print(f"    Saved to {args.output}")

    # Final summary
    print(f"\n{'='*60}")
    print(f"Summary — {args.output}")
    print(f"{'='*60}")
    print(f"  {'Key':45s} {'XS*lumi':>10} {'sum_abs_genw':>14} {'norm_w_mean':>14}")
    print("  " + "-"*87)
    for key in sorted(results):
        v = results[key]
        print(f"  {key:45s} {v['xs_lumi']:>10.1f} "
              f"{v['sum_abs_genw']:>14.4e} "
              f"{v['xs_lumi']/v['sum_abs_genw']:>14.6e}")


if __name__ == "__main__":
    main()
