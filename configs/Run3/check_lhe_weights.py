#!/usr/bin/env python3
"""
Generator weights of the Run 3 MC samples, from the first file of each dataset (datasets/<name>.json):

  * LHEPdfWeight: number of members and the branch title (it names the LHA IDs, i.e. the PDF set).
    pdf_uncertainty.py knows 103 (NNPDF31_nnlo_hessian_pdfas, 306000) and 101 (NNPDF31_nnlo_as_0118_nf_4,
    320900); the processor writes at most 103 members (custom_scale_factors.N_PDF_COLUMNS).
  * LHEScaleWeight / PSWeight: number of entries (sf_scaleweights expects 8 or 9, the PS weights 4).
  * mean |LHEWeight_originalXWGTUP|: the cross section per event only if the gridpack used event_norm =
    average; with event_norm = sum it is sigma / N (not usable as a cross section).
  * LHEScaleWeight grid: mean of each of the 9 weights and the fraction of events where it is 1 (the
    sample nominal); for the 4FS ttbb samples the nominal must sit at custom_scale_factors.
    TTBB_SAMPLE_NOMINAL_INDEX (5) and the analysis nominal (mu0, mu0) at TTBB_SCALE_REF_INDEX (4).
  * LHEReweightingWeight (with --reweighting, e.g. the tH ctcvcp samples): number of points and, per
    point, the mean weight = sigma(point) / sigma(generated). In Run 2 the tH samples were generated at
    kappa_t = -1 and the SM point (index 11 of 69) gave 0.074 (tHq) and 0.105 (tHW); if the 2024 samples
    are the same, the reweighting of Run2/th_reweighting.py is needed in Run 3 too.

    python3 check_lhe_weights.py                                # every MC sample of samples.DATASET_JSONS
    python3 check_lhe_weights.py TTX ttHTobb --reweighting      # selected samples, with the reweighting points
    python3 check_lhe_weights.py TTX --reweighting --datasets THQ_2024 THW_2024

Run where the files are readable (cluster). Only the first --events events of one file per dataset are read.
"""
import argparse
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KNOWN_PDF = {103: "NNPDF31_nnlo_hessian_pdfas (handled)", 101: "NNPDF31_nnlo_as_0118_nf_4 (handled)"}


def first_file(entry):
    files = entry.get("files") or []
    if not files:
        return None
    f = files[0]
    return f[len("file:"):] if f.startswith("file:") else f


def counts(arr):
    import awkward as ak
    n = np.unique(ak.to_numpy(ak.num(arr)))
    return n.tolist()


def check(path, n_events, reweighting):
    import uproot
    tree = uproot.open(path)["Events"]
    out = {}
    if "LHEWeight_originalXWGTUP" in tree:
        x = tree["LHEWeight_originalXWGTUP"].array(entry_stop=n_events, library="np")
        out["xwgtup"] = float(np.mean(np.abs(x)))
    if "LHEScaleWeight" in tree:
        sw = tree["LHEScaleWeight"].array(entry_stop=n_events)
        import awkward as ak
        if len(sw) and len(counts(sw)) == 1 and counts(sw)[0] == 9:
            w = ak.to_numpy(sw)
            out["scale_grid"] = [(float(np.mean(w[:, i])), float(np.mean(np.abs(w[:, i] - 1) < 5e-3))) for i in range(9)]
    for branch in ("LHEPdfWeight", "LHEScaleWeight", "PSWeight", "LHEReweightingWeight"):
        if branch not in tree:
            out[branch] = None
            continue
        arr = tree[branch].array(entry_stop=n_events)
        out[branch] = {"n": counts(arr), "title": tree[branch].title}
        if branch == "LHEReweightingWeight" and reweighting and len(arr):
            import awkward as ak
            if len(out[branch]["n"]) == 1 and out[branch]["n"][0] > 0:
                gen_sign = np.sign(tree["genWeight"].array(entry_stop=n_events, library="np"))
                w = ak.to_numpy(arr)
                out[branch]["mean"] = (w * gen_sign[:, None]).sum(axis=0) / gen_sign.sum()
    return out


def main():
    import samples
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="*", help="datasets/<name>.json names (default: all MC of samples.DATASET_JSONS)")
    ap.add_argument("--datasets", nargs="*", help="only these dataset names (e.g. THQ_2024)")
    ap.add_argument("--events", type=int, default=5000)
    ap.add_argument("--reweighting", action="store_true", help="mean weight per LHEReweightingWeight point")
    args = ap.parse_args()

    names = args.samples or [s for s in samples.DATASET_JSONS if s not in samples.DATA_SAMPLES]
    for name in names:
        path = os.path.join(HERE, "datasets", f"{name}.json")
        if not os.path.exists(path):
            print(f"{name}: datasets/{name}.json not built, skipped")
            continue
        for dataset, entry in json.load(open(path)).items():
            if args.datasets and dataset not in args.datasets:
                continue
            f = first_file(entry)
            print(f"\n=== {name} / {dataset}\n    {f}")
            if f is None:
                print("    no files")
                continue
            try:
                res = check(f, args.events, args.reweighting)
            except Exception as err:                      # noqa: BLE001  (report and go on)
                print(f"    cannot read: {err}")
                continue
            if "xwgtup" in res:
                print(f"    mean |LHEWeight_originalXWGTUP| = {res.pop('xwgtup'):.5g}  (sigma per event only with "
                      f"event_norm = average)")
            if "scale_grid" in res:
                grid = res.pop("scale_grid")
                print("    LHEScaleWeight [i] mean (fraction of events = 1):")
                print("        " + "  ".join(f"[{i}] {m:.4f} ({f:.2f})" for i, (m, f) in enumerate(grid)))
            for branch, info in res.items():
                if info is None:
                    print(f"    {branch:22s} absent")
                    continue
                note = ""
                if branch == "LHEPdfWeight":
                    n = info["n"]
                    note = KNOWN_PDF.get(n[0], "NOT HANDLED by pdf_uncertainty.py") if len(n) == 1 else "varying length"
                    if len(n) == 1 and n[0] > 103:
                        note += "; more than 103: the processor stops (N_PDF_COLUMNS)"
                print(f"    {branch:22s} entries {info['n']}  {note}\n        title: {info['title']}")
                if "mean" in info:
                    m = info["mean"]
                    print("        mean weight per point (sigma_i / sigma_generated):")
                    for i in range(0, len(m), 8):
                        print("        " + "  ".join(f"[{j:2d}] {m[j]:.4f}" for j in range(i, min(i + 8, len(m)))))


if __name__ == "__main__":
    main()
