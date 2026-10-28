#!/usr/bin/env python3
"""
MC b-tag efficiency maps of the Run 3 analysis, for the fixed-WP b-tag SFs (multi-WP method 1a,
custom_scale_factors.sf_btag_multiwp), from the btag_eff histogram of the btageff group:

    python3 run_pipeline.py inference --backend pbs --groups btageff
    python3 derive_btag_eff.py [-i outputs/nominal/btageff/output_all.coffea] [--plots plots/btageff]

Jets: JetGood of the analysis preselection, before the b-tag requirement (all jets whose b-tag
information the analysis uses). For each group of samples.BTAG_EFF_GROUPS (with the subsamples the
analysis takes from each sample), hadron flavour (b / c / light), pT and |eta| bin and working point
of analysis_years.BTAG_SF_WPS:

    eps_WP = sum of the weights of the jets tagged at WP (or tighter) / sum of the weights of all jets

Weights: genWeight x lumi x XS x pileup, normalized by the sum of genweights per dataset; they only set
the relative weight of the samples of a group. A bin of a group with fewer than --min-entries effective
jets takes the value of all groups together. eps is clipped to [0, 1] and made non-increasing from
the loosest to the tightest WP.

Output: a correctionlib file with one correction UParTAK4_eff_<group> per group, inputs
(working_point, flavor, abseta, pt), flow 'clamp' (pT above the last edge -> last bin). Then set
inference_common.BTAG_SFS[era]["eff"] to it.
"""
import argparse
import gzip
import os

import numpy as np
from coffea.util import load

import samples
from analysis_years import BTAG_SF_WPS

HERE = os.path.dirname(os.path.abspath(__file__))
HIST = "btag_eff"
CATEGORY = "presel"
TAGGER = "UParTAK4"          # pocket-coffea jet_scale_factors.btagSF[era].name
FLAVOURS = {0: "light", 4: "c", 5: "b"}


def dense(h):
    """(values, variances) over (flavor, pt, abseta, level) without flow bins, the pT overflow added
    to the last pT bin."""
    vals, var = h.values(flow=True), h.variances(flow=True)
    for i, ax in enumerate(h.axes):
        under = int(bool(getattr(ax.traits, "underflow", False)))
        n = len(ax)
        sl = [slice(None)] * vals.ndim
        if ax.name == "pt" and getattr(ax.traits, "overflow", False):
            last = [slice(None)] * vals.ndim
            last[i], over = under + n - 1, under + n
            ov = [slice(None)] * vals.ndim
            ov[i] = over
            vals[tuple(last)] += vals[tuple(ov)]
            var[tuple(last)] += var[tuple(ov)]
        sl[i] = slice(under, under + n)
        vals, var = vals[tuple(sl)], var[tuple(sl)]
    return vals, var


def group_histograms(output):
    """{group: (values, variances)} summed over the samples / subsamples / datasets of each group,
    and the axes of the histogram."""
    sums, axes = {}, None
    for key, datasets in output["variables"][HIST].items():
        sample, _, sub = key.partition("__")
        group = samples.btag_eff_group(sample, sub or None)
        if group is None:
            continue
        for dataset, h in datasets.items():
            h = h[{"cat": CATEGORY, "variation": "nominal"}]
            axes = axes or {ax.name: ax for ax in h.axes}
            vals, var = dense(h)
            if group in sums:
                sums[group] = (sums[group][0] + vals, sums[group][1] + var)
            else:
                sums[group] = (vals, var)
    if not sums:
        raise RuntimeError(f"No histogram '{HIST}' of a b-tag efficiency group in the input")
    return sums, axes


def efficiencies(vals, var):
    """eps[wp_index, flavor, pt, abseta] and the effective number of jets of the denominator."""
    total = vals.sum(axis=-1)
    total_var = var.sum(axis=-1)
    n_eff = np.where(total_var > 0, total ** 2 / np.where(total_var > 0, total_var, 1), 0)
    eps = []
    for k in range(1, len(BTAG_SF_WPS) + 1):
        tagged = vals[..., k:].sum(axis=-1)
        eps.append(np.where(total > 0, tagged / np.where(total > 0, total, 1), 0.0))
    eps = np.clip(np.array(eps), 0.0, 1.0)
    for k in range(1, len(eps)):
        eps[k] = np.minimum(eps[k], eps[k - 1])
    return eps, n_eff


def correction(name, eps, flavors, eta_edges, pt_edges):
    import correctionlib.schemav2 as cs

    def wp_node(k):
        return cs.Category(nodetype="category", input="flavor", content=[
            cs.CategoryItem(key=int(f), value=cs.MultiBinning(
                nodetype="multibinning", inputs=["abseta", "pt"], edges=[list(eta_edges), list(pt_edges)],
                content=[float(x) for x in eps[k, i].T.ravel()], flow="clamp"))
            for i, f in enumerate(flavors)])

    return cs.Correction(
        name=name, version=1,
        description=f"MC b-tag efficiency of {TAGGER} (fraction of jets tagged at the WP or tighter), "
                    f"derive_btag_eff.py",
        inputs=[cs.Variable(name="working_point", type="string", description="/".join(BTAG_SF_WPS)),
                cs.Variable(name="flavor", type="int", description="hadron flavour: 5=b, 4=c, 0=udsg"),
                cs.Variable(name="abseta", type="real"),
                cs.Variable(name="pt", type="real")],
        output=cs.Variable(name="efficiency", type="real"),
        data=cs.Category(nodetype="category", input="working_point", content=[
            cs.CategoryItem(key=wp, value=wp_node(k)) for k, wp in enumerate(BTAG_SF_WPS)]),
    )


def plot(eps, n_eff, group, flavors, eta_edges, pt_edges, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)
    centers = 0.5 * (pt_edges[1:] + pt_edges[:-1])
    for i, f in enumerate(flavors):
        fig, axs = plt.subplots(1, len(BTAG_SF_WPS), figsize=(5 * len(BTAG_SF_WPS), 4), sharey=False)
        for k, wp in enumerate(BTAG_SF_WPS):
            for j in range(len(eta_edges) - 1):
                axs[k].plot(centers, eps[k, i, :, j], marker="o",
                            label=rf"{eta_edges[j]:.1f} < |$\eta$| < {eta_edges[j + 1]:.1f}")
            axs[k].set_xscale("log")
            axs[k].set_xlabel(r"jet $p_T$ [GeV]")
            axs[k].set_ylabel(f"efficiency at {wp}")
            axs[k].set_title(f"{group}, {FLAVOURS.get(int(f), f)} jets, {wp}")
        axs[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"btag_eff_{group}_{FLAVOURS.get(int(f), f)}.png"), dpi=120)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", default=os.path.join(HERE, "outputs/nominal/btageff/output_all.coffea"))
    ap.add_argument("--era", default="2024")
    ap.add_argument("-o", "--output", default=None, help="default: btag_eff_<era>.json.gz in this directory")
    ap.add_argument("--min-entries", type=float, default=100.,
                    help="effective jets below which a group bin takes the value of all groups")
    ap.add_argument("--plots", default=None, help="directory for efficiency plots")
    args = ap.parse_args()
    out_path = args.output or os.path.join(HERE, f"btag_eff_{args.era}.json.gz")

    output = load(args.input)
    sums, axes = group_histograms(output)
    flavors = [int(lo) for lo in axes["flavor"].edges[:-1]]       # bins [0, 4), [4, 5), [5, 6)
    pt_edges, eta_edges = np.asarray(axes["pt"].edges), np.asarray(axes["abseta"].edges)
    if len(axes["level"]) != len(BTAG_SF_WPS) + 1:
        raise ValueError(f"Histogram has {len(axes['level'])} b-tag levels, expected {len(BTAG_SF_WPS) + 1} "
                         f"(analysis_years.BTAG_SF_WPS {BTAG_SF_WPS})")

    all_vals = sum(v for v, _ in sums.values())
    all_var = sum(v for _, v in sums.values())
    eps_all, neff_all = efficiencies(all_vals, all_var)
    low_all = int((neff_all < args.min_entries).sum())
    if low_all:
        print(f"  [warning] {low_all} bins have fewer than {args.min_entries:g} effective jets in all groups together")

    import correctionlib.schemav2 as cs
    corrections = []
    print(f"{'group':8s} {'flavour':7s} " + " ".join(f"{'eps_' + wp:>8s}" for wp in BTAG_SF_WPS)
          + f" {'jets (eff.)':>12s} {'bins from all':>14s}")
    for group in samples.BTAG_EFF_GROUPS:
        if group not in sums:
            print(f"  [warning] no jets for group {group}: it gets the map of all groups")
            eps, n_eff = eps_all, neff_all
        else:
            eps, n_eff = efficiencies(*sums[group])
        low = n_eff < args.min_entries
        eps = np.where(low[None], eps_all, eps)
        vals = sums.get(group, (all_vals, all_var))[0]
        for i, f in enumerate(flavors):
            # inclusive in pT / |eta|, for the summary
            tot = vals[i].sum()
            incl = [vals[i][..., k:].sum() / tot if tot > 0 else float("nan") for k in range(1, len(BTAG_SF_WPS) + 1)]
            print(f"{group:8s} {FLAVOURS.get(int(f), str(f)):7s} " + " ".join(f"{e:8.4f}" for e in incl)
                  + f" {n_eff[i].sum():12.0f} {int(low[i].sum()):8d}/{low[i].size}")
        corrections.append(correction(f"{TAGGER}_eff_{group}", eps, flavors, eta_edges, pt_edges))
        if args.plots:
            plot(eps, n_eff, group, flavors, eta_edges, pt_edges, args.plots)

    cset = cs.CorrectionSet(schema_version=2, corrections=corrections,
                            description=f"MC b-tag efficiency maps ({TAGGER}, {args.era}), derive_btag_eff.py")
    text = cset.model_dump_json(exclude_unset=True) if hasattr(cset, "model_dump_json") else cset.json(exclude_unset=True)
    with gzip.open(out_path, "wt") as fout:
        fout.write(text)
    print(f"\nWrote {out_path}: set inference_common.BTAG_SFS['{args.era}']['eff'] to it")


if __name__ == "__main__":
    main()
