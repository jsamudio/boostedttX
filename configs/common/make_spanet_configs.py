#!/usr/bin/env python3
"""
Write the configs of the SPANet input converters (../spanet_inputs) from the spanet_training.py of
an analysis (Run2 or Run3), run from that analysis directory:

    cd configs/Run2   (or Run3)
    python3 ../common/make_spanet_configs.py [-o spanet]

    spanet/features_parquet.yaml   coffea_to_parquet_script.py -c ...  (columns -> parquet)
    spanet/features_h5.yaml        parquet_to_h5.py -c ...             (parquet -> SPANet h5)
    spanet/event_info_inputs.yaml  the INPUTS section of the SPANet event-info file

It also checks that the SPANet inputs can be fed by the inference (spanet_training.INFERENCE_FEATURES)
and prints the spanet_models.py entry of a model trained on these inputs.
"""
import argparse
import os
import sys

import yaml

sys.path.insert(0, os.getcwd())          # the analysis directory, with its spanet_training.py
import spanet_training as st  # noqa: E402

# SPANet normalization of the inputs (event-info file); every other feature: "normalize"
LOG_NORMALIZE = {"pt"}
NO_NORMALIZATION = {"btagB", "bbtag", "is_electron", "xbb_wp_L", "xbb_wp_M", "xbb_wp_T"}
SEQUENTIAL_INPUTS = {"Jet"}


def parquet_config():
    features, pad = {}, {}
    for collection, spec in st.EXPORT.items():
        features[collection] = dict(spec["fields"])
        names = set(spec["fields"])
        # the converter zips the collections as Momentum4D: give them a mass (and eta for the MET)
        if "pt" in names and "phi" in names:
            missing = {"mass": 0} if not names & {"mass", "m", "energy", "E"} else {}
            if "eta" not in names:
                missing["eta"] = 0
            if missing:
                pad[collection] = missing
    return {
        "features": {"common": features, "by_sample": {}},
        "features_pad": {"common": pad, "by_sample": {}},
        "awkward_collections": [c for c, spec in st.EXPORT.items() if spec.get("jagged")],
        "matched_collections": dict(st.MATCHED_COLLECTIONS),
    }


def h5_config():
    classes = list(dict.fromkeys(label for label, _ in st.CLASSES.values()))
    return {
        "input": {name: list(features) for name, (_, features) in st.SPANET_INPUTS.items()},
        "collection": {name: collection for name, (collection, _) in st.SPANET_INPUTS.items()},
        "particles": list(st.PARTICLES),
        "mapping_sample": {key: label for key, (label, _) in st.CLASSES.items()},
        "mapping_encoding": {"signal": {key: index for key, (_, index) in st.CLASSES.items()}},
        "classification": {"EVENT": classes + ["signal"]},
        "frac_train": st.FRAC_TRAIN,
    }


def event_info_inputs():
    """INPUTS section of the SPANet event-info file, from SPANET_INPUTS."""
    def norm(feature):
        if feature in LOG_NORMALIZE:
            return "log_normalize"
        return "none" if feature in NO_NORMALIZATION else "normalize"
    out = {"SEQUENTIAL": {}, "GLOBAL": {}}
    for name, (_, features) in st.SPANET_INPUTS.items():
        kind = "SEQUENTIAL" if name in SEQUENTIAL_INPUTS else "GLOBAL"
        out[kind][name] = {f: norm(f) for f in features if f != "MASK"}
    return {"INPUTS": out}


def check_inference():
    """Warnings for SPANet inputs the inference cannot provide."""
    warnings = []
    for name, (_, features) in st.SPANET_INPUTS.items():
        provided = st.INFERENCE_FEATURES.get(name)
        if provided is None:
            warnings.append(f"{name}: input not fed by the inference")
        elif name in st.INFERENCE_FIXED and list(features) != list(provided):
            warnings.append(f"{name}: {features} != {provided} used by the inference")
        elif set(features) - set(provided) - {"MASK"}:
            warnings.append(f"{name}: {sorted(set(features) - set(provided) - {'MASK'})} not provided by the inference")
    return warnings


def check_classes():
    """Each CLASSES key must match exactly one (sub)sample name; indices must be 0..N-1."""
    names = []
    for sample, opts in st.SAMPLES.items():
        subs = opts.get("subsamples")
        if subs == "ttbar":
            names += [f"{sample}__{s}" for s in ("tt+LF", "tt+C", "tt+B")]
        elif subs:
            names += [f"{sample}__{s}" for s in subs]
        else:
            names.append(sample)
    problems = []
    for name in names:
        keys = [k for k in st.CLASSES if k in name]
        if len(keys) != 1:
            problems.append(f"{name}: matches {keys or 'no'} CLASSES key(s)")
    indices = sorted({index for _, index in st.CLASSES.values()})
    if indices != list(range(len(indices))):
        problems.append(f"class indices {indices} are not 0..{len(indices) - 1}: SPANet would get "
                        f"{max(indices) + 1} classes, some of them empty")
    return names, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default=os.path.join(os.getcwd(), "spanet"))
    args = ap.parse_args()

    st.validate()
    names, problems = check_classes()
    if problems:
        sys.exit("CLASSES problems:\n  " + "\n  ".join(problems))

    os.makedirs(args.outdir, exist_ok=True)
    header = "# Generated by make_spanet_configs.py from spanet_training.py: edit that file and rerun.\n"
    outputs = {"features_parquet.yaml": parquet_config(), "features_h5.yaml": h5_config(),
               "event_info_inputs.yaml": event_info_inputs()}
    for name, cfg in outputs.items():
        path = os.path.join(args.outdir, name)
        with open(path, "w") as f:
            f.write(header)
            yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=None)
        print(f"wrote {path}")

    for w in check_inference():
        print(f"WARNING (inference): {w}")

    jet = [f for f in st.SPANET_INPUTS["Jet"][1] if f != "MASK"]
    fatjet = list(st.SPANET_INPUTS["FatJet"][1])
    event = list(st.SPANET_INPUTS["Event"][1]) if "Event" in st.SPANET_INPUTS else []
    print("\n(Sub)samples -> one parquet each:", ", ".join(names))
    print("spanet_models.py entry for a model trained on these inputs:")
    print(f'    "<name>": {{"path": "<model>.onnx",\n'
          f'               "jet_features": {jet},\n'
          f'               "fatjet_features": {fatjet},\n'
          f'               "event_features": {event},\n'
          f'               "count_clean_jets": {st.COUNT_CLEAN_JETS},\n'
          f'               "clean_jets": {st.CLEAN_JETS}}},')


if __name__ == "__main__":
    main()
