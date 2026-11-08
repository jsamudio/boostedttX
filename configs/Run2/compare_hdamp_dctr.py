#!/usr/bin/env python3
"""
hdamp: dedicated 4FS ttbb hdamp samples vs DCTR reweighting of the nominal ttbb (Run 2, 2018).

Question: does the DCTR hdamp reweighting (classifier on the gen top kinematics, trained on ttbar)
reproduce the hdamp up / down effect on the 4FS ttbb in this selection? If so, it can stand in for
the hdamp samples that do not exist for the Run 3 ttbb.

Inputs (.coffea, analysis selection, any era):
  --nominal   nominal TTbb_* run with PC_DCTR=1 (columns dctr_hdamp_up/down + output['dctr_sums'])
  --hdamp     TTbb_*_hdampUP / _hdampDOWN run (group 'variations')

Normalization (both keep the total cross section of the nominal sample):
  hdamp samples: lumi * xsec(nominal) / sum sign(genWeight) of the variation sample (outputToPkl.get_norm)
  DCTR:          event weight * dctr_hdamp_<var> * S_nom / S_<var>, S = sum sign(genWeight) * w over ALL
                 generated events of the dataset (dctr_sums 'all'), so only the selection acceptance moves.
Event weight = datacard weight (datacard_shapes.getZhbbWeight), selection = makeDatacard.cuts.

Output: per subsample (default tt+B) the rate ratio to nominal (total and per Z/H pT channel) for
hdamp sample vs DCTR, and per observable the ratio to nominal of both, with a shape chi2 of
(sample - DCTR) using the stat. uncertainties (approximate: nominal and DCTR share the events), and
a 'flat' chi2 of the sample against a pure rate change (shape ratio = 1), which is the argument for
treating hdamp as a rate-only uncertainty (as in the previous analysis).
Plots: <outdir>/<subsample>_<observable>.png/.pdf

    python3 compare_hdamp_dctr.py --nominal outputs_dctr/nominal/bkg/output_all.coffea \
        --hdamp outputs_hdamp/nominal/variations/output_all.coffea [--btag-norm btag_norm_hdamp.json]
"""
import argparse
import gc
import os

import numpy as np
import pandas as pd

import outvars
from samples import TTBAR_SPLIT
from spanet_models import ACTIVE_MODELS, score_columns

TTBB = ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu"]
VARS = {"up": "hdampUP", "down": "hdampDOWN"}
CHANNELS = {"Zhpt1": (200, 300), "Zhpt2": (300, 450), "Zhpt3": (450, np.inf)}

OBSERVABLES = {          # column: (bins, label)
    "ZH_pt": (np.array([200, 250, 300, 375, 450, 525, 600, 700, 850, 1200]), r"Z/H candidate $p_T$ [GeV]"),
    "ZH_M": (np.linspace(50, 200, 16), "Z/H candidate mass [GeV]"),
    "ZH_bbvLscore": (np.linspace(0.9, 1.0, 11), "Xbb score"),
    "n_ak4jets": (np.arange(4.5, 12.5, 1.0), "N AK4 jets"),
    "n_b_inZH": (np.arange(-0.5, 4.5, 1.0), "N b jets in Z/H"),
    "n_q_outZH": (np.arange(-0.5, 6.5, 1.0), "N light jets outside Z/H"),
    "MET_pt": (np.array([20, 50, 80, 120, 170, 250, 400, 800]), r"$p_T^{miss}$ [GeV]"),
}
for _m in ACTIVE_MODELS:
    OBSERVABLES[f"disc_{_m}"] = (np.linspace(0, 1, 11), f"disc {_m}")
    OBSERVABLES[f"ttbb_score_{_m}"] = (np.linspace(0, 1, 11), f"ttbb score {_m}")
    OBSERVABLES[f"ttlf_score_{_m}"] = (np.linspace(0, 1, 11), f"ttlf score {_m}")


def read(files, samples, btag_norm, dctr=False, era=None):
    """DataFrame of the selected events of the given (sub)samples; with dctr=True also the DCTR
    normalization factors S_nom / S_var per dataset (columns dctr_norm_up/down). era: keep only that era
    (e.g. the 2018 part of a multi-era production output)."""
    from coffea.util import load
    from outputToPkl import OutputToPkl
    from datacard_shapes import getZhbbWeight
    from makeDatacard import cuts

    reader = OutputToPkl([], btag_norm)
    variables = (outvars.common_vars + outvars.weight_vars
                 + [c for m in ACTIVE_MODELS for c in score_columns(m)])
    frames = []
    for path in files:
        print(f"Reading {path}")
        filein = load(path)
        genweight_df = pd.DataFrame.from_dict(filein["sum_signOf_genweights"], orient="index")
        for sample in samples:
            for df in reader.read_samples(filein, genweight_df, [sample], variables):
                if era is not None and df["era"].iloc[0] != era:
                    continue
                df["subsample"] = sample.split("__")[1]
                if dctr:
                    base = sample.split("__")[0]
                    sums = filein.get("dctr_sums", {}).get(base, {})
                    if len(sums) != 1:
                        raise KeyError(f"expected one dataset in dctr_sums[{base}] (PC_DCTR=1 run, one era), "
                                       f"got {list(sums)}")
                    s = next(iter(sums.values()))["all"]
                    for v in VARS:
                        df[f"dctr_norm_{v}"] = s["nominal"] / s[f"dctr_hdamp_{v}"]
                frames.append(df)
        del filein
        gc.collect()
    if not frames:
        raise ValueError(f"no events of {samples} in {files}")
    df = pd.concat(frames, ignore_index=True)
    df = df[cuts(df, str(df["year"].iloc[0]))].copy()
    df["w"] = getZhbbWeight(df).to_numpy(dtype=float)
    for m in ACTIVE_MODELS:
        s, b, lf = (df[f"{c}_score_{m}"] for c in ("sig", "ttbb", "ttlf"))
        den = (s + b + lf).to_numpy(dtype=float)
        df[f"disc_{m}"] = np.where(den > 0, s / np.where(den > 0, den, 1.0), 0.0)
    return df


def hist(x, w, bins):
    x = np.clip(np.asarray(x, dtype=float), bins[0], bins[-1] - 1e-9)       # under/overflow into the edges
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w ** 2)
    return h, np.sqrt(e2)


def ratio(h, e, hn, en):
    r = np.divide(h, hn, out=np.zeros_like(h), where=hn > 0)
    re = np.sqrt(np.divide(e, hn, out=np.zeros_like(h), where=hn > 0) ** 2
                 + (r * np.divide(en, hn, out=np.zeros_like(h), where=hn > 0)) ** 2)
    return r, re


def yields(df, w, mask=None):
    m = np.ones(len(df), bool) if mask is None else mask
    return w[m].sum(), np.sqrt((w[m] ** 2).sum())


def chi2_pvalue(c2, ndf):
    """p-value of chi2 for ndf degrees of freedom (scipy if available), as a short string."""
    if ndf <= 0:
        return "-"
    try:
        from scipy.stats import chi2
        return f"{chi2.sf(c2, ndf):.2f}"
    except ImportError:
        return "n/a"


def compare(nom, var, sub, outdir, plot=True):
    nom, var = nom[nom["subsample"] == sub], {v: d[d["subsample"] == sub] for v, d in var.items()}
    wn = nom["w"].to_numpy()
    wd = {v: wn * nom[f"dctr_hdamp_{v}"].to_numpy() * nom[f"dctr_norm_{v}"].to_numpy() for v in VARS}
    ws = {v: var[v]["w"].to_numpy() for v in VARS}

    print(f"\n=== {sub}: {len(nom)} nominal events, " + ", ".join(f"{VARS[v]} {len(var[v])}" for v in VARS))
    for v in VARS:
        dw = nom[f"dctr_hdamp_{v}"].to_numpy()
        print(f"  DCTR {v}: event weight mean {dw.mean():.4f} min {dw.min():.3f} max {dw.max():.3f}, "
              f"gen-level norm S_nom/S_var = {nom[f'dctr_norm_{v}'].iloc[0]:.4f}")
    print(f"  {'region':8} {'nominal':>10}" + "".join(f" {VARS[v] + ' / nom':>20} {'DCTR ' + v + ' / nom':>20}"
                                                      for v in VARS))
    regions = {"all": (None, None)}
    regions.update({ch: lim for ch, lim in CHANNELS.items()})
    for name, (lo, hi) in regions.items():
        mn = None if lo is None else ((nom["ZH_pt"] > lo) & (nom["ZH_pt"] <= hi)).to_numpy()
        yn, en = yields(nom, wn, mn)
        line = f"  {name:8} {yn:10.2f}"
        for v in VARS:
            mv = None if lo is None else ((var[v]["ZH_pt"] > lo) & (var[v]["ZH_pt"] <= hi)).to_numpy()
            ys, es = yields(var[v], ws[v], mv)
            yd, ed = yields(nom, wd[v], mn)
            rs, rse = ratio(np.array([ys]), np.array([es]), np.array([yn]), np.array([en]))
            # DCTR / nominal: same events, so the stat. uncertainty is that of the weight difference
            dd = np.sqrt((((wd[v] - wn)[mn] if mn is not None else (wd[v] - wn)) ** 2).sum())
            line += f" {rs[0]:12.4f} +- {rse[0]:.4f} {yd / yn:12.4f} +- {dd / yn:.4f}"
        print(line)

    print(f"  shape chi2 / ndf (sample vs DCTR, both normalized to the nominal yield;"
          f" 'flat' = sample vs a pure rate change, i.e. sample/nominal at unit area against 1):")
    for obs, (bins, label) in OBSERVABLES.items():
        if obs not in nom:
            continue
        hn, en = hist(nom[obs], wn, bins)
        curves, chi = {}, []
        for v in VARS:
            hs, es = hist(var[v][obs], ws[v], bins)
            hd, _ = hist(nom[obs], wd[v], bins)
            hdiff, ediff = hist(nom[obs], wd[v] - wn, bins)
            # shape: scale each to the nominal yield
            ks, kd = hn.sum() / hs.sum(), hn.sum() / hd.sum()
            rs, rse = ratio(hs * ks, es * ks, hn, en)
            rd = np.divide(hd * kd, hn, out=np.zeros_like(hn), where=hn > 0)
            rde = np.divide(ediff, hn, out=np.zeros_like(hn), where=hn > 0)
            ok = (hn > 0) & (hs > 0)
            c2 = (((rs - rd) ** 2) / np.maximum(rse ** 2 + rde ** 2, 1e-12))[ok].sum()
            # pure-rate hypothesis: the shape ratio is flat at 1 (sample and nominal stat., independent);
            # one dof less for the normalization
            c2f = (((rs - 1.0) ** 2) / np.maximum(rse ** 2, 1e-12))[ok].sum()
            chi.append(f"{v} {c2:.1f}/{ok.sum() - 1} (flat {c2f:.1f}, p={chi2_pvalue(c2f, ok.sum() - 1)})")
            curves[v] = (rs, rse, rd, rde)
        print(f"    {obs:28} " + "  ".join(chi))
        if plot:
            draw(bins, label, curves, f"{outdir}/{sub.replace('+', '')}_{obs}")


def draw(bins, label, curves, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    x = 0.5 * (bins[1:] + bins[:-1])
    colors = {"up": "#d62728", "down": "#1f77b4"}
    for v, (rs, rse, rd, rde) in curves.items():
        ax.errorbar(x, rs, yerr=rse, fmt="o", ms=4, color=colors[v], label=f"{VARS[v]} sample / nominal")
        ax.stairs(rd, bins, baseline=None, color=colors[v], lw=1.5, label=f"DCTR {v} / nominal")
        ax.fill_between(bins, np.append(rd - rde, rd[-1] - rde[-1]), np.append(rd + rde, rd[-1] + rde[-1]),
                        step="post", color=colors[v], alpha=0.2, lw=0)
    ax.axhline(1, color="gray", lw=0.8)
    ax.set_xlabel(label)
    ax.set_ylabel("ratio to nominal (shape)")
    vals = np.concatenate([np.concatenate([rs[rs > 0], rd[rd > 0]]) for rs, _, rd, _ in curves.values()])
    span = max(0.1, np.abs(vals - 1).max() if len(vals) else 0.1)
    ax.set_ylim(1 - 1.15 * min(span, 0.5), 1 + 1.6 * min(span, 0.5))       # headroom for the legend
    ax.legend(fontsize=7, ncol=2, loc="upper center")
    ax.set_title("4FS tt+bb, 2018: hdamp samples vs DCTR", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)

DC_MASS_BINS = [50, 80, 105, 145, 200]          # makeDatacard.mass_bins


def datacard_bins(nom, var, sub, edges_file, model, outdir, plot=True):
    """hdamp sample / nominal per datacard bin (NN bin x mass bin, flattened as NN-major: bin = i_nn * 4 + i_mass),
    per pT channel, with the NN edges of the cards (nn_edges_disc_<model>.json). Prints the inclusive rate,
    the flat chi2 (shape vs pure rate) in the datacard binning, and draws the per-bin ratio with the
    nominal stat. (band) and the hdamp-sample stat. (error bars), as in the previous AN's hdamp figures."""
    import json
    with open(edges_file) as f:
        all_edges = json.load(f)
    year = str(nom["year"].iloc[0])
    edges = all_edges.get(year, next(iter(all_edges.values())))
    col = f"disc_{model}"
    nom = nom[nom["subsample"] == sub]
    var = {v: d[d["subsample"] == sub] for v, d in var.items()}
    print(f"\n  datacard binning ({col}, edges {edges_file}, {year}): hdamp sample / nominal")
    res = {}
    for ch, (lo, hi) in CHANNELS.items():
        nn_edges = np.asarray(edges[ch], dtype=float)

        def h2(df):
            m = ((df["ZH_pt"] > lo) & (df["ZH_pt"] <= hi)).to_numpy()
            x, mass, w = (df[c].to_numpy(dtype=float)[m] for c in (col, "ZH_M", "w"))
            h, _, _ = np.histogram2d(x, mass, bins=[nn_edges, DC_MASS_BINS], weights=w)
            e2, _, _ = np.histogram2d(x, mass, bins=[nn_edges, DC_MASS_BINS], weights=w ** 2)
            return h.ravel(), np.sqrt(e2.ravel())

        hn, en = h2(nom)
        res[ch] = {"nominal": (hn, en)}
        line = f"    {ch}: {len(hn)} bins"
        for v in VARS:
            hv, ev = h2(var[v])
            rate = hv.sum() / hn.sum()
            ok = (hn > 0) & (hv > 0)
            # shape: scale the variation to the nominal total, ratio per bin, stat. of both
            r = np.divide(hv / rate, hn, out=np.ones_like(hn), where=hn > 0)
            re = np.sqrt(np.divide(ev / rate, hn, out=np.zeros_like(hn), where=hn > 0) ** 2
                         + (r * np.divide(en, hn, out=np.zeros_like(hn), where=hn > 0)) ** 2)
            c2 = (((r - 1) / np.maximum(re, 1e-12)) ** 2)[ok].sum()
            line += f"   {VARS[v]}: rate {rate:.4f}, flat chi2 {c2:.1f}/{ok.sum() - 1} (p={chi2_pvalue(c2, ok.sum() - 1)})"
            res[ch][v] = (hv, ev, rate)
        print(line)
    # one lnN = one common rate in every bin of every channel: flat chi2 over all bins with the global rate
    hn_all = np.concatenate([res[ch]["nominal"][0] for ch in res])
    en_all = np.concatenate([res[ch]["nominal"][1] for ch in res])
    line = f"    all channels ({len(hn_all)} bins, one common rate = lnN):"
    for v in VARS:
        hv = np.concatenate([res[ch][v][0] for ch in res])
        ev = np.concatenate([res[ch][v][1] for ch in res])
        rate = hv.sum() / hn_all.sum()
        ok = (hn_all > 0) & (hv > 0)
        r = np.divide(hv / rate, hn_all, out=np.ones_like(hn_all), where=hn_all > 0)
        re = np.sqrt(np.divide(ev / rate, hn_all, out=np.zeros_like(hn_all), where=hn_all > 0) ** 2
                     + (r * np.divide(en_all, hn_all, out=np.zeros_like(hn_all), where=hn_all > 0)) ** 2)
        c2 = (((r - 1) / np.maximum(re, 1e-12)) ** 2)[ok].sum()
        line += f"   {VARS[v]}: rate {rate:.4f}, flat chi2 {c2:.1f}/{ok.sum() - 1} (p={chi2_pvalue(c2, ok.sum() - 1)})"
    print(line)
    if plot:
        draw_datacard_bins(res, sub, f"{outdir}/{sub.replace('+', '')}_datacard_bins")


def draw_datacard_bins(res, sub, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(res), 1, figsize=(8, 3.2 * len(res)), squeeze=False)
    colors = {"up": "#1f77b4", "down": "#ff7f0e"}
    for ax, (ch, r) in zip(axes[:, 0], res.items()):
        hn, en = r["nominal"]
        n = len(hn)
        edges = np.arange(n + 1)
        rel = np.divide(en, hn, out=np.zeros_like(hn), where=hn > 0)
        ax.fill_between(edges, np.append(1 - rel, 1 - rel[-1]), np.append(1 + rel, 1 + rel[-1]), step="post",
                        facecolor="none", edgecolor="#d62728", hatch="xxxx", lw=0, label="nominal stat. unc.")
        for k, v in enumerate(VARS):
            hv, ev, rate = r[v]
            ratio_ = np.divide(hv, hn, out=np.full_like(hn, np.nan), where=hn > 0)
            err = np.divide(ev, hn, out=np.zeros_like(hn), where=hn > 0)
            ax.errorbar(np.arange(n) + 0.35 + 0.3 * k, ratio_, yerr=err, fmt="o", ms=3, color=colors[v],
                        label=f"{VARS[v]} / nominal (sample stat.)")
            ax.axhline(rate, color=colors[v], ls="--", lw=1.2, label=f"{VARS[v]} rate {rate:.3f}")
        ax.axhline(1, color="gray", lw=0.8)
        ax.set_xlim(0, n)
        ax.set_ylim(0.0, 2.0)
        ax.set_ylabel("Variation / Nom")
        ax.set_title(f"tt+B hdamp, {ch} (bin = NN bin x 4 + mass bin)", fontsize=9)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, ncol=3, loc="upper left")
    axes[-1, 0].set_xlabel("Bin number")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    plt.close(fig)



if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nominal", nargs="+", required=True, help="nominal TTbb .coffea with PC_DCTR=1")
    ap.add_argument("--hdamp", nargs="+", required=True, help="TTbb hdampUP/DOWN .coffea")
    ap.add_argument("--btag-norm", default=None,
                    help="btag_norm.json with entries for the nominal AND hdamp TTbb samples (btag_norm.py on both "
                         "outputs); default: b-tag SF not normalized (same for both, cancels to first order)")
    ap.add_argument("--subsamples", nargs="+", default=["tt+B"], choices=list(TTBAR_SPLIT))
    ap.add_argument("--outdir", default="plots/hdamp_dctr")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--edges", default=None,
                    help="NN edges of the cards (nn_edges_disc_<model>.json from makeDatacard.py): also compare "
                         "per datacard bin (pT channel x NN bin x mass bin)")
    ap.add_argument("--model", default=None,
                    help="SPANet model of --edges (default: taken from the file name nn_edges_disc_<model>.json)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    nom = read(args.nominal, [f"{s}__{c}" for s in TTBB for c in args.subsamples], args.btag_norm, dctr=True)
    var = {v: read(args.hdamp, [f"{s}_{tag}__{c}" for s in TTBB for c in args.subsamples], args.btag_norm)
           for v, tag in VARS.items()}
    for sub in args.subsamples:
        compare(nom, var, sub, args.outdir, plot=not args.no_plots)
        if args.edges:
            model = args.model or os.path.basename(args.edges).removeprefix("nn_edges_disc_").removesuffix(".json")
            datacard_bins(nom, var, sub, args.edges, model, args.outdir, plot=not args.no_plots)
