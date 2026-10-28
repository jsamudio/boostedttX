#!/usr/bin/env python3
"""
pT spectra of the EFT tt+bb sample for the SM and for one Wilson coefficient at several values
(all others 0), with the EFT/SM ratio:

  * gen bb pT, all generated events  - from the accumulator in the .coffea output (eft_genbb.py)
  * gen bb pT, selected events       - from the dumped columns (events_gen_bb_pt)
  * reco Z/H candidate pT, selected  - from the dumped columns (events_ZH_pt)

The selected events are those of the EFT dump (btag_mask category); --datacard-cuts adds the
datacard selection (>= 5 AK4 jets, 2 b jets outside the candidate, Xbb >= the analysis working
point of --era, MET >= 20 GeV, 50 < mass < 200 GeV).

The EFT weight of an event at a WC point c (c_0 = 1 for the SM) is
    w(c) = sum_k EFTfitCoefficients_k * c[EFTfitCoefficientIndex1_k] * c[EFTfitCoefficientIndex2_k]

With --xsec (and --coffea, which holds the SM weight sum of all generated events) the spectra are
expected events: w(c) * xsec * lumi / (sum of the SM weights of all generated events).

    python3 plot_eft_pt.py --wc cbwre --values -2 -1 1 2 \\
        --coffea outputs_eft/nominal/EFT/output_all.coffea \\
        --dumps /cms/data/store/user/jsamudio/eftchunks/2018/ttbbSMEFT_2018/btag_mask/nominal \\
        --datacard-cuts --era 2018 --xsec 19.261 --lumi 59.8 -o plots/eft_pt
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from eft_genbb import GEN_BB_PT_BINS

# variable: (dump column, bins, x label)
VARIABLES = {
    "gen_bb_pt": ("events_gen_bb_pt", GEN_BB_PT_BINS, r"$p_{T}^{b\bar{b},\ gen}$ [GeV]"),
    "ZH_pt":     ("events_ZH_pt",     np.arange(200, 1050, 50), r"$p_{T}^{Z/H\ cand.}$ [GeV]"),
}
EFT_DUMP_COLUMNS = ["events_EFTfitCoefficients", "events_EFTfitCoefficientIndex1",
                    "events_EFTfitCoefficientIndex2", "events_WCnames"]
CUT_COLUMNS = ["events_n_ak4jets", "events_n_b_outZH", "events_ZH_bbvLscore", "events_MET_pt", "events_ZH_M"]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#e600e6", "#8c564b", "#b2182b"]


def decode_WCnames(WCnames):
    """WC names from the integer encoding of the EFT samples (as in EFTScaling.ipynb)."""
    out, name = [], ""
    for val in WCnames:
        four = int(val).to_bytes(4, "big").decode()
        if "-" in four:
            name += four.split("-")[-1]
        else:
            if name:
                out.append(name)
            name = four.lstrip("\x00")
    out.append(name)
    return out


def point_weights(index1, index2, names, wc, value):
    """Vector over the coefficients: c[index1] * c[index2] at the point wc = value."""
    if wc not in names:
        raise SystemExit(f"{wc} is not one of the sample's WCs: {names[1:]}")
    c = np.array([1.0] + [value if n == wc else 0.0 for n in names[1:]])
    return c[np.asarray(index1)] * c[np.asarray(index2)]


def binned_sums(values, coefs, bins):
    """(n_bins, n_coefficients) sums of the coefficient vectors, values clipped into the last bin."""
    sums = np.zeros((len(bins) - 1, coefs.shape[1]))
    ok = values >= bins[0]
    idx = np.clip(np.digitize(values[ok], bins) - 1, 0, len(bins) - 2)
    np.add.at(sums, idx, coefs[ok])
    return sums


# ---------------------------------------------------------------- inputs
def read_accumulator(paths):
    """Gen bb pT coefficient sums over all generated events:
    (sums, index1, index2, names, coefficient sums of all generated events or None)."""
    from coffea.util import load
    total, structure, all_events = None, None, None
    for path in paths:
        out = load(path)
        if "eft_genbb" not in out:
            raise SystemExit(f"No 'eft_genbb' in {path}: produced before eft_genbb.py was added?")
        for dataset, per_structure in out["eft_genbb"].items():
            for key, sums in per_structure.items():
                if structure is not None and key != structure:
                    raise SystemExit(f"{dataset}: different EFT coefficient structure than the other datasets")
                structure = key
                total = sums if total is None else total + sums
            for key, tot in out.get("eft_genbb_total", {}).get(dataset, {}).items():
                all_events = tot if all_events is None else all_events + tot
    meta = json.loads(structure)
    if all_events is None:
        print("  [WARNING] no 'eft_genbb_total' in the .coffea output (produced before it was added): "
              "--xsec normalization not possible")
    return total, meta["index1"], meta["index2"], ["SM"] + decode_WCnames(meta["wcnames"]), all_events


def read_dumps(paths, datacard_cuts, era):
    """Coefficient sums per variable from the dumped columns: ({variable: sums}, index1, index2, names)."""
    import awkward as ak
    import pyarrow.dataset as ds
    wp = None
    if datacard_cuts:
        from analysis_years import bbvl_wp
        wp = bbvl_wp(era)
    sums, meta = {}, None
    for path in paths:
        dataset = ds.dataset(path, format="parquet")
        available = set(dataset.schema.names)
        variables = {v: spec for v, spec in VARIABLES.items() if spec[0] in available}
        missing = [spec[0] for v, spec in VARIABLES.items() if v not in variables]
        if missing:
            print(f"  [{path}] no {missing}: those spectra are skipped")
        columns = EFT_DUMP_COLUMNS + [s[0] for s in variables.values()] + (CUT_COLUMNS if datacard_cuts else [])
        for batch in dataset.to_batches(columns=columns):
            ch = ak.from_arrow(batch)
            if len(ch) == 0:
                continue
            if meta is None:
                meta = (ak.to_numpy(ch["events_EFTfitCoefficientIndex1"][0]).astype(int),
                        ak.to_numpy(ch["events_EFTfitCoefficientIndex2"][0]).astype(int),
                        ["SM"] + decode_WCnames(ak.to_list(ch["events_WCnames"][0])))
            coefs = ak.to_numpy(ch["events_EFTfitCoefficients"]).astype(np.float64)
            keep = np.ones(len(ch), dtype=bool)
            if datacard_cuts:
                col = {c: ak.to_numpy(ch[c]) for c in CUT_COLUMNS}
                keep = ((col["events_n_ak4jets"] >= 5) & (col["events_n_b_outZH"] == 2)
                        & (col["events_ZH_bbvLscore"] >= wp) & (col["events_MET_pt"] >= 20)
                        & (col["events_ZH_M"] >= 50) & (col["events_ZH_M"] <= 200))
            for v, (column, bins, _) in variables.items():
                s = binned_sums(ak.to_numpy(ch[column])[keep], coefs[keep], bins)
                sums[v] = s if v not in sums else sums[v] + s
    if meta is None:
        raise SystemExit("No events in the dumps")
    return sums, meta[0], meta[1], meta[2]


# ---------------------------------------------------------------- plot
def plot(sums, index1, index2, names, wc, values, bins, xlabel, label, out, normalize, scale=1.0):
    plt.style.use(hep.style.CMS)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    sums = sums * scale
    sm = sums @ point_weights(index1, index2, names, wc, 0.0)
    norm = (lambda h: h / h.sum()) if normalize else (lambda h: h)
    hep.histplot(norm(sm), bins, color="black", lw=2, label="SM", ax=ax)
    for value, color in zip(values, COLORS):
        spec = sums @ point_weights(index1, index2, names, wc, value)
        hep.histplot(norm(spec), bins, color=color, lw=2, label=f"{wc} = {value:g}", ax=ax)
        ratio = np.divide(norm(spec), norm(sm), out=np.full_like(sm, np.nan), where=norm(sm) != 0)
        rax.stairs(ratio, bins, color=color, lw=2)
    rax.axhline(1, color="black", lw=1)
    rax.set_ylabel("EFT / SM")
    rax.set_xlabel(xlabel)
    ax.set_ylabel("Normalized" if normalize else ("Events / bin" if scale != 1.0 else "Weighted events (a.u.)"))
    ax.set_yscale("log")
    rax.set_ylim(0,4)
    ax.set_xlim(bins[0], bins[-1])
    ax.legend(title=r"$t\bar{t}+b\bar{b}$ EFT sample, " + label, fontsize="small", title_fontsize="small")
    hep.cms.label("Preliminary", data=False, rlabel="(13 TeV)", ax=ax)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(os.path.splitext(out)[0] + ".png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wc", required=True, help="Wilson coefficient to vary (others at 0)")
    ap.add_argument("--values", nargs="+", type=float, required=True)
    ap.add_argument("--coffea", nargs="*", default=[], help=".coffea output(s) with eft_genbb (all generated events)")
    ap.add_argument("--dumps", nargs="*", default=[], help="parquet dump directories of the EFT sample")
    ap.add_argument("--datacard-cuts", action="store_true", help="apply the datacard selection to the dumps")
    ap.add_argument("--era", default="2018", help="era of the Xbb working point for --datacard-cuts")
    ap.add_argument("--normalize", action="store_true", help="compare shapes (each spectrum to unit area)")
    ap.add_argument("--xsec", type=float, default=None, help="SM cross section [pb] -> expected events (needs --coffea)")
    ap.add_argument("--lumi", type=float, default=59.8, help="luminosity [fb^-1] for --xsec")
    ap.add_argument("-o", "--outdir", default="plots/eft_pt")
    args = ap.parse_args()
    if not args.coffea and not args.dumps:
        raise SystemExit("Give --coffea and/or --dumps")
    if args.xsec is not None and not args.coffea:
        raise SystemExit("--xsec needs --coffea (the SM weight sum of all generated events)")

    tag = f"{args.wc}{'_norm' if args.normalize else ''}"
    scale = 1.0
    if args.coffea:
        sums, i1, i2, names, all_events = read_accumulator(args.coffea)
        if args.xsec is not None and all_events is not None:
            k00 = [k for k, (a, b) in enumerate(zip(i1, i2)) if a == 0 and b == 0][0]
            scale = args.xsec * args.lumi * 1000 / all_events[k00]
            print(f"  expected events: xsec {args.xsec} pb x {args.lumi} fb^-1 / SM weight sum "
                  f"{all_events[k00]:.4g} -> scale {scale:.4g}")
        plot(sums, i1, i2, names, args.wc, args.values, GEN_BB_PT_BINS, VARIABLES["gen_bb_pt"][2],
             "all generated events", os.path.join(args.outdir, f"gen_bb_pt_inclusive_{tag}.pdf"),
             args.normalize, scale)
    if args.dumps:
        sums, i1, i2, names = read_dumps(args.dumps, args.datacard_cuts, args.era)
        label = "datacard selection" if args.datacard_cuts else "b-tag category"
        for v, s in sums.items():
            _, bins, xlabel = VARIABLES[v]
            plot(s, i1, i2, names, args.wc, args.values, bins, xlabel, label,
                 os.path.join(args.outdir, f"{v}_selected_{tag}.pdf"), args.normalize, scale)


if __name__ == "__main__":
    main()