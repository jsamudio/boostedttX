#!/usr/bin/env python3
"""
Training vs validation loss per epoch of a SPANet run, from its TensorBoard event files: the classic
overtraining check on the loss curves (the validation loss reaches a minimum and rises while the
training loss keeps falling).

Curves (per epoch):
    classification   training 'loss/classification/<head>' (per step, averaged over the epoch) and
                     'validation_classification_loss' (per epoch): the same cross entropy (with the
                     class balancing weights), the checkpoint / early-stopping metric of options_run3.json
    total            training 'loss/total_loss' / number of loss terms (the objective is their mean;
                     the logged value is their sum) and 'validation_loss' (mean of the terms)
    assignment       training 'loss/<particle>/assignment_loss' (no validation counterpart is logged)
Dropout is on in training and off in validation, so the validation loss can sit below the training
loss; overtraining shows as the validation curve turning up (or the gap growing) while the training
curve still falls, not as the size of the gap. The best checkpoint (best_checkpoint.json) and the
last epoch are marked.

    python3 plot_training_losses.py -l <SPANet>/run3_refactor/version_0 [-o plots/training]

Needs the tensorboard package (the SPANet environment has it).
"""
import argparse
import glob
import json
import os
import re

import numpy as np


def load_scalars(log_dir):
    from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
    files = sorted(glob.glob(os.path.join(log_dir, "events.out.tfevents.*")))
    if not files:
        raise SystemExit(f"no TensorBoard event files in {log_dir}")
    scalars = {}
    for path in files:
        acc = EventAccumulator(path, size_guidance={"scalars": 0})
        acc.Reload()
        for tag in acc.Tags()["scalars"]:
            ev = acc.Scalars(tag)
            steps, vals = scalars.get(tag, ([], []))
            scalars[tag] = (steps + [e.step for e in ev], vals + [e.value for e in ev])
    return {t: (np.asarray(s), np.asarray(v, dtype=float)) for t, (s, v) in scalars.items()}


def to_epochs(steps, values, epoch_steps, epoch_values):
    """Epoch of each logged step (Lightning's 'epoch' scalar), then the mean per epoch."""
    order = np.argsort(epoch_steps)
    es, ev = epoch_steps[order], epoch_values[order]
    idx = np.clip(np.searchsorted(es, steps, side="right") - 1, 0, len(es) - 1)
    epochs = np.rint(ev[idx]).astype(int)
    uniq = np.unique(epochs)
    return uniq, np.array([values[epochs == e].mean() for e in uniq])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-l", "--log-dir", required=True, help="SPANet run directory (version_N)")
    ap.add_argument("-o", "--outdir", default="plots/training")
    args = ap.parse_args()

    sc = load_scalars(args.log_dir)
    if "epoch" not in sc:
        raise SystemExit(f"no 'epoch' scalar in {args.log_dir}; tags: {sorted(sc)}")
    e_steps, e_vals = sc["epoch"]
    per_epoch = lambda tag: to_epochs(*sc[tag], e_steps, e_vals)

    best_epoch, best_info = None, ""
    record = os.path.join(args.log_dir, "best_checkpoint.json")
    if os.path.exists(record):
        best = json.load(open(record))
        m = re.search(r"epoch=(\d+)", os.path.basename(best.get("path", "")))
        best_epoch = int(m.group(1)) if m else None
        best_info = f"best checkpoint: epoch {best_epoch} ({best.get('metric')} = {best.get('score'):.4f})"

    opts = {}
    if os.path.exists(os.path.join(args.log_dir, "options.json")):
        opts = json.load(open(os.path.join(args.log_dir, "options.json")))

    panels = []
    cls_tags = [t for t in sc if t.startswith("loss/classification/")]
    if cls_tags and "validation_classification_loss" in sc:
        train = [per_epoch(t) for t in cls_tags]                     # one head normally; sum the heads
        epochs = train[0][0]
        tr = np.sum([np.interp(epochs, e, v) for e, v in train], axis=0) * float(opts.get("classification_loss_scale", 1.0))
        panels.append(("classification loss (cross entropy)", (epochs, tr), per_epoch("validation_classification_loss")))
    if "loss/total_loss" in sc and "validation_loss" in sc:
        n_terms = None
        # terms of the objective: one assignment (+ detection) per particle, symmetric KL, classification
        part = {t.split("/")[1] for t in sc if re.match(r"loss/[^/]+/assignment_loss$", t)}
        n_terms = len(part) * (2 if any(re.match(r"loss/[^/]+/detection_loss$", t) for t in sc) else 1) \
            + int("loss/symmetric_loss" in sc) + len(cls_tags) \
            + len([t for t in sc if t.startswith("loss/regression/")])
        e, v = per_epoch("loss/total_loss")
        panels.append((f"total loss (mean of {n_terms} terms)", (e, v / max(n_terms, 1)), per_epoch("validation_loss")))
    assign = [t for t in sc if re.match(r"loss/[^/]+/assignment_loss$", t)]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(args.outdir, exist_ok=True)
    n = len(panels) + int(bool(assign))
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 4.2), squeeze=False)
    print(f"{args.log_dir}: {best_info}")
    for ax, (title, (te, tv), (ve, vv)) in zip(axes[0], panels):
        ax.plot(te, tv, color="#1f77b4", lw=1.6, label="training (epoch mean, dropout on)")
        ax.plot(ve, vv, color="#2ca02c", lw=1.6, marker="o", ms=2.5, label="validation")
        i = int(np.argmin(vv))
        ax.axvline(ve[i], color="#2ca02c", ls=":", lw=1, label=f"validation minimum (epoch {ve[i]})")
        if best_epoch is not None:
            ax.axvline(best_epoch, color="black", ls="--", lw=1, label=f"best checkpoint (epoch {best_epoch})")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("epoch")
        ax.set_ylabel("loss")
        ax.legend(fontsize=7)
        last = vv[-1]
        rise = (last - vv[i]) / vv[i] if vv[i] else np.nan
        tr_min = np.interp(ve[i], te, tv)
        print(f"  {title}: validation min {vv[i]:.4f} at epoch {ve[i]}, last {last:.4f} (epoch {ve[-1]}, "
              f"{100 * rise:+.1f}% vs min); training at the minimum {tr_min:.4f}, last {tv[-1]:.4f} "
              f"({100 * (tv[-1] - tr_min) / tr_min:+.1f}%)")
    if assign:
        ax = axes[0][-1]
        for t in assign:
            e, v = per_epoch(t)
            ax.plot(e, v, lw=1.4, label=f"training {t.split('/')[1]}")
        if best_epoch is not None:
            ax.axvline(best_epoch, color="black", ls="--", lw=1)
        ax.set_title("assignment loss (training only)", fontsize=10)
        ax.set_xlabel("epoch")
        ax.legend(fontsize=7)
    fig.suptitle(f"SPANet training: {os.path.normpath(args.log_dir)}", fontsize=10, x=0.02, ha="left")
    fig.tight_layout()
    stem = os.path.join(args.outdir, "training_losses")
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    print(f"Wrote {stem}.png/.pdf")


if __name__ == "__main__":
    main()
