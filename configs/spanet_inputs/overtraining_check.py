#!/usr/bin/env python3
"""
SPANet overtraining check: classification scores of a trained model on its training events against
events it never saw (the test h5 of parquet_to_h5.py), per true class.

Three sets, all from the files of the run (log directory options.json):
    train       training file, range [0, split)           (split = dataset_limit x train_validation_split)
    validation  training file, range [split, dataset_limit) (selects the checkpoint and stops the training)
    test        --test (the h5 converter's test file; must not be the run's validation_file)
The ranges are those SPANet used (contiguous, dataset_randomization must be 0).

Per true class (CLASSIFICATIONS/EVENT/<target>) and per score (each class probability and the analysis
discriminant disc = (p0 + p1) / sum, signal = classes --signal-classes): unit-area weighted
distributions (WEIGHTS/weight), train (filled) vs test (points) and validation (line), ratio to test;
weighted chi2/ndf over --bins bins (stat. of both) and the unweighted two-sample KS p-value, for
train vs test and for validation vs test. Validation vs test is the reference of two sets the model
was not trained on: with millions of events a KS p-value is small for tiny differences, so overtraining
shows as train vs test clearly worse than validation vs test, not as a small p-value alone.

    cd SPANet   (environment with torch and spanet)
    python3 ../configs/spanet_inputs/overtraining_check.py -l spanet_output/version_0 \
        -t <dir>/run3_test_<N>.h5 [--max-events 300000] [--gpu] [-o plots/overtraining]

Panels: each true class, the merged signal (--signal-classes, e.g. ttH + ttZ, same disc shape: twice the
statistics) and the merged background; in a merged group each class has the same total weight. Threshold table: weighted fraction of each group with
--threshold-score (disc) >= --thresholds, train / validation / test, with uncertainties, train/test
ratio and pull: a training-specific gain in high-score efficiency (what would bias the templates)
shows as positive pulls for train and not for validation.

Outputs: <outdir>/overtraining_<score>.png/.pdf (one panel per group), overtraining_summary.csv
(chi2, KS) and overtraining_thresholds.csv; both tables are also printed.
"""
import argparse
import csv
import json
import os

import h5py
import numpy as np

CLASS_NAMES = ["tthbb", "ttzbb", "ttbb", "ttlf"]      # spanet_training.CLASSES order (index 0..3)


def evaluate(model, path, lo, hi, max_events, gpu):
    """Softmax scores (n, n_classes), true class (n,) and weights (n,) of the events [lo, hi) of an h5
    file, as SPANet indexes them (fractions of the file; contiguous). At most max_events events."""
    import torch
    from spanet.evaluation import evaluate_on_test_dataset

    with h5py.File(path, "r") as f:
        n_file = len(f["WEIGHTS/weight"])
        lower, upper = int(round(lo * n_file)), int(round(hi * n_file))
        if max_events and upper - lower > max_events:
            upper = lower + max_events
        weights = f["WEIGHTS/weight"][lower:upper].astype(float)
    model.testing_dataset = model.dataset(
        data_file=path, event_info=model.options.event_info_file,
        limit_index=(lower / n_file, upper / n_file), vector_limit=model.options.limit_to_num_jets)
    ds = model.testing_dataset
    if len(ds) != upper - lower:
        raise RuntimeError(f"{path}: SPANet selected {len(ds)} events, expected {upper - lower}")
    keys = [k for k in ds.classifications if k.endswith("signal")] or list(ds.classifications)
    labels = ds.classifications[keys[0]].numpy().astype(int)
    if gpu:
        model = model.cuda()
    with torch.no_grad():
        evaluation = evaluate_on_test_dataset(model, progress=None)
    cls = evaluation.classifications
    key = [k for k in cls if k.endswith("signal")] or list(cls)
    scores = np.asarray(cls[key[0]])
    print(f"  {os.path.basename(path)} [{lower}, {upper}): {len(scores)} events")
    return scores, labels, weights


def hist(x, w, bins):
    h, _ = np.histogram(np.clip(x, bins[0], bins[-1] - 1e-9), bins, weights=w)
    e2, _ = np.histogram(np.clip(x, bins[0], bins[-1] - 1e-9), bins, weights=w ** 2)
    s = h.sum()
    return (h / s, np.sqrt(e2) / s) if s > 0 else (h, np.sqrt(e2))


def group_weights(y, w, classes):
    """Mask and weights of the events of a group of true classes. A merged group gives each class the same
    total weight (each class normalized to 1): the check is about the score shape of each class, and with
    the physics weights one class can dominate (e.g. tt+LF in the background)."""
    m = np.isin(y, classes)
    if len(classes) == 1:
        return m, w[m]
    wg = w.astype(float).copy()
    for c in classes:
        mc = y == c
        tot = wg[mc].sum()
        if tot != 0:
            wg[mc] = wg[mc] / tot
    return m, wg[m]


def efficiency(passed, w):
    """Weighted fraction passing and its uncertainty (independent events: the variance of
    sum_pass w / sum w is (sum_pass w^2 (1 - e)^2 + sum_fail w^2 e^2) / (sum w)^2)."""
    tot = w.sum()
    if tot <= 0:
        return np.nan, np.nan
    e = w[passed].sum() / tot
    var = ((w[passed] ** 2).sum() * (1 - e) ** 2 + (w[~passed] ** 2).sum() * e ** 2) / tot ** 2
    return float(e), float(np.sqrt(var))


def chi2(a, b):
    (ha, ea), (hb, eb) = a, b
    den = ea ** 2 + eb ** 2
    ok = den > 0
    return float((((ha - hb) ** 2)[ok] / den[ok]).sum()), int(ok.sum()) - 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-l", "--log-dir", required=True, help="SPANet run directory (options.json, checkpoints)")
    ap.add_argument("-t", "--test", required=True, help="held-out h5 (parquet_to_h5.py test file)")
    ap.add_argument("-c", "--checkpoint", default=None, help="default: the run's best checkpoint")
    ap.add_argument("--max-events", type=int, default=300000, help="per set (0 = all)")
    ap.add_argument("--signal-classes", nargs="+", type=int, default=[0, 1], help="class indices of the disc numerator")
    ap.add_argument("--classes", nargs="+", default=CLASS_NAMES, help="class names by index")
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--ratio-range", nargs=2, type=float, default=[0.7, 1.3])
    ap.add_argument("--threshold-score", default="disc", help="score of the threshold table")
    ap.add_argument("--thresholds", nargs="+", type=float, default=[0.5, 0.7, 0.8, 0.9, 0.95])
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("-o", "--outdir", default="plots/overtraining")
    args = ap.parse_args()

    from scipy.stats import ks_2samp
    from spanet.evaluation import load_model

    opts = json.load(open(os.path.join(args.log_dir, "options.json")))
    if opts.get("validation_file"):
        print(f"[warning] the run used validation_file = {opts['validation_file']}: if that is --test, the test set "
              "also selected the checkpoint (biased check); train with validation from the training file")
        if os.path.abspath(opts["validation_file"]) == os.path.abspath(args.test):
            print("[warning] --test IS the run's validation file")
    if int(opts.get("dataset_randomization", 0)) != 0:
        raise SystemExit("dataset_randomization != 0: the training range is not contiguous, not supported")
    split = float(opts.get("dataset_limit", 1.0)) * float(opts.get("train_validation_split", 0.95))
    limit = float(opts.get("dataset_limit", 1.0))
    train_file = opts["training_file"]

    model = load_model(args.log_dir, testing_file=args.test, checkpoint=args.checkpoint, cuda=args.gpu)
    sets = {"train": evaluate(model, train_file, 0.0, split, args.max_events, args.gpu)}
    if not opts.get("validation_file"):
        sets["validation"] = evaluate(model, train_file, split, limit, args.max_events, args.gpu)
    sets["test"] = evaluate(model, args.test, 0.0, 1.0, args.max_events, args.gpu)

    n_cls = sets["test"][0].shape[1]
    names = (list(args.classes) + [f"class{i}" for i in range(n_cls)])[:n_cls]

    def score_values(s):
        out = {f"p_{names[i]}": s[:, i] for i in range(n_cls)}
        tot = s.sum(axis=1)
        out["disc"] = np.divide(s[:, args.signal_classes].sum(axis=1), tot, out=np.zeros(len(s)), where=tot > 0)
        return out
    values = {k: (score_values(s), y, w) for k, (s, y, w) in sets.items()}
    bins = np.linspace(0, 1, args.bins + 1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(args.outdir, exist_ok=True)
    # panels: each true class, then the merged signal (--signal-classes) and the merged background
    others = [i for i in range(n_cls) if i not in args.signal_classes]
    groups = [(f"true {names[c]}", [c]) for c in range(n_cls)]
    groups += [("signal (" + " + ".join(names[i] for i in args.signal_classes) + ")", list(args.signal_classes)),
               ("background (" + " + ".join(names[i] for i in others) + ")", others)]
    ncol = 3
    nrow = (len(groups) + ncol - 1) // ncol
    rows = []
    for score in values["test"][0]:
        fig = plt.figure(figsize=(4.9 * ncol, 4.75 * nrow))
        ratios = []
        for r in range(nrow):
            ratios += [3, 1] + ([0.55] if r < nrow - 1 else [])
        gs = fig.add_gridspec(len(ratios), ncol, height_ratios=ratios, hspace=0.06, wspace=0.24)
        for g, (gname, gcls) in enumerate(groups):
            row = 3 * (g // ncol)
            ax = fig.add_subplot(gs[row, g % ncol])
            rax = fig.add_subplot(gs[row + 1, g % ncol], sharex=ax)
            ax.tick_params(labelbottom=False)
            h = {}
            for k, (v, y, w) in values.items():
                m, wg = group_weights(y, w, gcls)
                h[k] = hist(v[score][m], wg, bins)
            x = 0.5 * (bins[1:] + bins[:-1])
            ax.stairs(h["train"][0], bins, fill=True, alpha=0.35, color="#1f77b4", label="train")
            if "validation" in h:
                ax.stairs(h["validation"][0], bins, color="#2ca02c", lw=1.4, label="validation")
            ax.errorbar(x, h["test"][0], yerr=h["test"][1], fmt="ko", ms=3, label="test")
            ht = h["test"][0]
            ok = ht > 0
            for k, color in (("train", "#1f77b4"), ("validation", "#2ca02c")):
                if k in h:
                    r = np.divide(h[k][0], ht, out=np.full_like(ht, np.nan), where=ok)
                    rax.stairs(r, bins, baseline=None, color=color, lw=1.4)
            rax.axhline(1, color="gray", lw=0.8)
            rax.set_ylim(*args.ratio_range)
            rax.set_ylabel("/ test", fontsize=8)
            text = []
            for k in ("train", "validation"):
                if k not in h:
                    continue
                c2, ndf = chi2(h[k], h["test"])
                (vk, yk, _), (vt, yt, _) = values[k], values["test"]
                mk, mt = np.isin(yk, gcls), np.isin(yt, gcls)
                p = ks_2samp(vk[score][mk], vt[score][mt]).pvalue
                text.append(f"{k} vs test: $\\chi^2$/ndf {c2:.0f}/{ndf}, KS p {p:.2g}")
                rows.append({"score": score, "true_class": gname.replace("true ", ""), "comparison": f"{k}_vs_test",
                             "chi2": round(c2, 2), "ndf": ndf, "ks_p": p,
                             "n_" + k: int(mk.sum()), "n_test": int(mt.sum())})
            ax.set_title(gname, fontsize=9)
            ax.text(0.02, 0.97, "\n".join(text), transform=ax.transAxes, va="top", fontsize=7,
                    bbox={"facecolor": "white", "alpha": 0.7, "lw": 0})
            ax.set_yscale("log")
            pos = ht[ht > 0]
            if len(pos):
                ax.set_ylim(0.3 * pos.min(), 30 * pos.max())
            ax.set_ylabel("unit area", fontsize=8)
            rax.set_xlabel(score)
            if g == 0:
                handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right", ncol=3, fontsize=8, frameon=False)
        fig.suptitle(f"SPANet overtraining check: {score} ({os.path.basename(os.path.normpath(args.log_dir))})",
                     fontsize=10, x=0.02, ha="left")
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(args.outdir, f"overtraining_{score}.{ext}"))
        plt.close(fig)

    # fraction of events above thresholds of --threshold-score, train / validation / test, per group
    trows = []
    for gname, gcls in groups:
        for t in args.thresholds:
            row = {"group": gname.replace("true ", ""), "score": args.threshold_score, "threshold": t}
            for k, (v, y, w) in values.items():
                m, wg = group_weights(y, w, gcls)
                e, de = efficiency(v[args.threshold_score][m] >= t, wg)
                row[f"eff_{k}"], row[f"err_{k}"] = e, de
            et, dt = row["eff_test"], row["err_test"]
            for k in ("train", "validation"):
                if f"eff_{k}" not in row:
                    continue
                ek, dk = row[f"eff_{k}"], row[f"err_{k}"]
                row[f"{k}_over_test"] = ek / et if et > 0 else np.nan
                row[f"pull_{k}"] = (ek - et) / np.hypot(dk, dt) if np.hypot(dk, dt) > 0 else np.nan
            trows.append(row)
    has_val = "validation" in values
    print(f"\nFraction with {args.threshold_score} >= threshold (weighted; pull = (eff - eff_test) / combined error)")
    print(f"{'group':>32} {'cut':>5} {'train':>15} " + (f"{'validation':>15} " if has_val else "")
          + f"{'test':>15} {'train/test':>10} {'pull':>6}" + (f" {'val pull':>8}" if has_val else ""))
    for r in trows:
        fmt = lambda k: f"{r['eff_' + k]:.4f}+-{r['err_' + k]:.4f}"
        print(f"{r['group']:>32} {r['threshold']:5.2f} {fmt('train'):>15} " + (f"{fmt('validation'):>15} " if has_val else "")
              + f"{fmt('test'):>15} {r['train_over_test']:10.3f} {r['pull_train']:6.2f}"
              + (f" {r['pull_validation']:8.2f}" if has_val else ""))
    tpath = os.path.join(args.outdir, "overtraining_thresholds.csv")
    with open(tpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(trows[0]))
        writer.writeheader()
        writer.writerows(trows)
    print(f"Wrote {tpath}")

    print(f"\n{'score':>12} {'true':>26} {'comparison':>20} {'chi2/ndf':>12} {'KS p':>9}")
    for r in rows:
        print(f"{r['score']:>12} {r['true_class']:>26} {r['comparison']:>20} {r['chi2']:>7.1f}/{r['ndf']:<4} {r['ks_p']:9.2g}")
    path = os.path.join(args.outdir, "overtraining_summary.csv")
    with open(path, "w", newline="") as f:
        fields = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {args.outdir}/overtraining_<score>.png/.pdf and {path}")


if __name__ == "__main__":
    main()
