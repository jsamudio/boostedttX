#!/usr/bin/env python3
"""
List the corrections and correction sets used by a pocket-coffea run, per era, from the files
pocket-coffea writes into every output directory:
    parameters_dump.yaml   all resolved parameters (cvmfs files with their version folder, keys)
    config.json            the workflow options (the analysis's own correction files)

    python3 ../common/list_corrections.py outputs/nominal/signal --eras 2016_PreVFP 2016_PostVFP 2017 2018
    python3 ../common/list_corrections.py outputs/nominal/signal --eras 2024 > corrections_2024.txt

Any parameter keyed by era is reduced to the requested era; parameters that exist only for other
eras are left out.
"""
import argparse
import json
import os

import yaml

ALL_ERAS = ["2016_PreVFP", "2016_PostVFP", "2017", "2018",
            "2022_preEE", "2022_postEE", "2023_preBPix", "2023_postBPix", "2024", "2025"]

# Sections of parameters_dump.yaml that hold corrections
SECTIONS = {
    "Jet energy corrections / resolution": ["jets_calibration", "jet_types"],
    "JEC reapplied to data / MC": ["jets_calibration", "apply_jec_Data"],
    "Jet scale factors (b tag, PU ID, veto maps, jet ID)": ["jet_scale_factors"],
    "b-tagging working points": ["btagging", "working_point"],
    "Lepton scale factors and momentum corrections": ["lepton_scale_factors"],
    "Pileup": ["pileupJSONfiles"],
    "MET filters": ["event_flags"],
    "MET filters (data only)": ["event_flags_data"],
    "Luminosity": ["lumi", "picobarns"],
}
# Workflow options that point to the analysis's own correction files
OWN_FILES = ["sf_dir"]


def for_era(node, era):
    """Reduce every dict keyed by eras to the entry of `era` (None if that era is absent)."""
    if isinstance(node, dict):
        if any(k in ALL_ERAS for k in node):
            return for_era(node[era], era) if era in node else None
        reduced = {k: for_era(v, era) for k, v in node.items()}
        return {k: v for k, v in reduced.items() if v is not None} or None
    return node


def get(node, path):
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", help="pocket-coffea output directory (with parameters_dump.yaml)")
    ap.add_argument("--eras", nargs="+", required=True)
    args = ap.parse_args()

    params = yaml.safe_load(open(os.path.join(args.outdir, "parameters_dump.yaml")))
    config_path = os.path.join(args.outdir, "config.json")
    options = {}
    if os.path.exists(config_path):
        options = json.load(open(config_path)).get("workflow", {}).get("workflow_options", {})

    for era in args.eras:
        print("#" * 80 + f"\n# {era}\n" + "#" * 80)
        for title, path in SECTIONS.items():
            node = for_era(get(params, path), era)
            print(f"\n== {title}")
            print(yaml.safe_dump(node, default_flow_style=False, width=200, sort_keys=False).rstrip()
                  if node is not None else "  (none)")
        print("\n== Analysis's own correction files (workflow options)")
        for key in OWN_FILES:
            print(f"  {key}: {options.get(key, '(not set)')}")
        print("  (sf_dir holds ele_trig_sf_<era>.json.gz, muon_trig_sf_<era>.json.gz and the Xbb SF files)\n")


if __name__ == "__main__":
    main()
