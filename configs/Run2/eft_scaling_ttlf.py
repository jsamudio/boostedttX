#!/usr/bin/env python3
"""
EFT scaling (Combine interferenceModel JSON) for the SPANet categorizations of makeDatacard.py.
Two-step (--ttlf-cut <cut> [--ttlf-cr]):
    SR: ttlf score <  cut, NN axis s / (s + ttbb), edges nn_edges_sigvsttbb_<model>_ttlf<cut>.json
    CR: ttlf score >= cut, one bin of the disc s / (s + ttbb + ttlf)
Argmax (--argmax [--argmax-lf]), categories by the highest score (datacard_shapes.spanet_category):
    ZH:  highest = sig,  NN axis the disc,       edges nn_edges_disc_<model>_argmaxZH.json
    TTB: highest = ttbb, NN axis the ttbb score, edges nn_edges_ttbb_score_<model>_argmaxTTB.json
    LF:  highest = ttlf, one bin of the disc
Same method and defaults as EFTScaling.ipynb (Cells 1 and 2): the 2018 SMEFT dumps, the datacard
selection with the 2018 Xbb LP, per-bin EFT/SM ratios with the low-statistics fallback
(channel-level, then process-level ratio) and the diagonal mass-peak bins kept. The EFT samples
exist for 2018 only: one scaling file per year of the edges file, each with that year's edges.

Usage (after makeDatacard.py --ttlf-cut 0.2 --ttlf-cr, from configs/Run2):
    python3 eft_scaling_ttlf.py --cut 0.2
    python3 eft_scaling_ttlf.py --cut 0.2 --no-cr          # SR cards only
    python3 eft_scaling_ttlf.py --argmax [--no-cr]         # after makeDatacard.py --argmax [--argmax-lf]

Output: combine_scaling_ttlf<cut>_<year>.json with SR and CR entries, or combine_scaling_argmax_<year>.json
with ZH, TTB and LF entries (--no-cr: no CR / LF). Channel names follow --channel-format, default
"{region}_{channel}" (e.g. SR_Zhpt1, TTB_Zhpt2), which matches
    combineCards.py SR=<SR card> CR=<CR card>      or      combineCards.py ZH=<ZH> TTB=<TTB> [LF=<LF>]
Change it if the cards are combined with other labels.
"""
import argparse
import concurrent.futures
import json

import numpy as np

# ---------------- Analysis configuration: keep in sync with the datacard (as EFTScaling.ipynb) ----------------
GEN_PT_BINS = [0, 200, 300, 450]       # signal split on genZHpt: [0,200), [200,300), [300,450), >=450
MASS_BINS = [50, 80, 105, 145, 200]
BBVL_WP = 0.9172                       # datacard: ZH_bbvLscore >= LP of 2018 (the EFT samples are 2018)
# Reco channels, right-inclusive like pd.cut in the datacard: (200,300], (300,450], (450,inf)
CHANNELS = {"Zhpt1": (200, 300), "Zhpt2": (300, 450), "Zhpt3": (450, np.inf)}
MERGE_LAST_MASS_BIN = {"Zhpt1"}        # datacard merges the last two mass bins in Zhpt1 only (isblind=True)

BASE = "/cms/data/store/user/jsamudio/eftchunks/2018"
PROCESSES = {
    "ttH": {"paths": [f"{BASE}/ttHSMEFT_2018/btag_mask/nominal"], "gen_split": True},
    "ttZ": {"paths": [f"{BASE}/ttZSMEFT_2018/btag_mask/nominal"], "gen_split": True},
    "tt_B": {"paths": [f"{BASE}/ttbbSMEFT_2018/btag_mask/nominal"], "gen_split": False},
}
# Wilson coefficients and range of the nominal scaling files (EFTScaling.ipynb, Cell 2)
DEFAULT_WCS = ["chq1", "chq3", "cht", "chtbre", "ctwre", "ctbre", "cthre", "cbwre", "cqq8", "cqtqb1re",
               "ctb1", "cqb1", "ctb8", "cqtqb8re", "cqq1", "cqt1", "cqt8"]
WC_RANGE = (-40, 40)

DIAGONAL = {"Zhpt1": 1, "Zhpt2": 2, "Zhpt3": 3}   # reco channel -> matching gen-pT bin index
BOSON_MASS = {"ttH": 125.0, "ttZ": 91.19}
CATEGORIES = ("ZH", "TTB", "LF")            # order of the scores, as datacard_shapes.spanet_category


def score_columns(model):
    return [f"events_{s}_score_{model}" for s in ("sig", "ttbb", "ttlf")]


def base_columns(model):
    return ["events_EFTfitCoefficients", "events_EFTfitCoefficientIndex1",
            "events_EFTfitCoefficientIndex2", "events_WCnames",
            "events_n_ak4jets", "events_n_b_outZH", "events_ZH_bbvLscore",
            "events_MET_pt", "events_ZH_M", "events_ZH_pt", *score_columns(model)]


# ---------------- Helpers (as EFTScaling.ipynb) ----------------
def decode_WCnames(WCnames):
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


def bin_index(values, edges):
    """Bin like np.histogram: left-inclusive, last bin also includes the upper edge."""
    edges = np.asarray(edges, dtype=float)
    idx = np.searchsorted(edges, values, side="right") - 1
    idx[values == edges[-1]] = len(edges) - 2
    return idx


def fill(acc, nn, mass, coefs, nn_edges, mass_edges):
    i, j = bin_index(nn, nn_edges), bin_index(mass, mass_edges)
    ok = (i >= 0) & (i < len(nn_edges) - 1) & (j >= 0) & (j < len(mass_edges) - 1)
    np.add.at(acc, (i[ok], j[ok]), coefs[ok])


def gen_masks(label, gen_split, genpt, n):
    if not gen_split:
        return {label: np.ones(n, dtype=bool)}
    edges = GEN_PT_BINS + [np.inf]
    return {f"{label}{i}": (genpt >= edges[i]) & (genpt < edges[i + 1]) for i in range(len(GEN_PT_BINS))}


def unroll(w, channel):
    """Same layout as the datacard: (nn, mass) flattened row-major, mass fastest."""
    if channel in MERGE_LAST_MASS_BIN:
        w = w.copy()
        w[:, -2, :] += w[:, -1, :]
        w = w[:, :-1, :]
    return w.reshape(-1, w.shape[-1])


def mass_edges_of(channel):
    mass = list(MASS_BINS)
    return mass[:-2] + [mass[-1]] if channel in MERGE_LAST_MASS_BIN else mass


def diagonal_mask(channel, process, nn_edges):
    """Bins of a diagonal entry (gen-pT bin == reco channel) in the boson's mass-peak bin."""
    mass = mass_edges_of(channel)
    n_mass = len(mass) - 1
    n = (len(nn_edges) - 1) * n_mass
    i = DIAGONAL.get(channel)
    if i is None or not process[-1:].isdigit() or int(process[-1]) != i:
        return np.zeros(n, dtype=bool)
    m0 = BOSON_MASS.get(process.rstrip("0123456789"))
    if m0 is None:
        return np.zeros(n, dtype=bool)
    j_peak = int(np.searchsorted(mass, m0, side="right") - 1)
    return np.arange(n) % n_mass == j_peak


def unreliable(n_raw, sm, sm2, min_events, max_rel_err):
    """True where the EFT-sample statistics are too poor for a per-bin ratio."""
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_err = np.where(sm > 0, np.sqrt(sm2) / sm, np.inf)
    return (n_raw < min_events) | (sm <= 0) | (rel_err > max_rel_err)


# ---------------- Reading the EFT samples ----------------
def region_scores(col, model, cut):
    """{region: (event mask, NN value)}. cut = ttlf cut: SR = ttlf < cut in s/(s+ttbb), CR = ttlf >= cut
    in the disc. cut = 'argmax': ZH / TTB / LF = highest score sig / ttbb / ttlf (ties to the lower index),
    in the disc / the ttbb score / the disc."""
    s, b, lf = (col(c).astype(np.float64) for c in score_columns(model))
    den_sr, den_cr = s + b, s + b + lf
    sr = np.where(den_sr > 0, s / np.where(den_sr > 0, den_sr, 1.0), 0.0)
    cr = np.where(den_cr > 0, s / np.where(den_cr > 0, den_cr, 1.0), 0.0)
    if cut == "argmax":
        cat = np.argmax(np.column_stack([s, b, lf]), axis=1)
        return {"ZH": (cat == 0, cr), "TTB": (cat == 1, b), "LF": (cat == 2, cr)}
    return {"SR": (lf < cut, sr), "CR": (lf >= cut, cr)}


def process_path(label, cfg, path, model, cut, edges):
    """Accumulate binned EFT coefficients (+ raw count and SM weight^2) for one directory.
    edges: {region: {year: {channel: nn edges}}}. Returns ({region: {year: {(channel, process): array}}}, meta)."""
    import awkward as ak
    import pyarrow.dataset as ds
    print(f"[{label}] reading {path}")
    columns = base_columns(model) + (["events_genZHpt"] if cfg["gen_split"] else [])
    acc = {r: {y: {} for y in e} for r, e in edges.items()}
    meta, k00 = None, None
    for batch in ds.dataset(path, format="parquet").to_batches(columns=columns):
        chunk = ak.from_arrow(batch)
        if len(chunk) == 0:
            continue
        if meta is None:
            idx1 = ak.to_numpy(chunk["events_EFTfitCoefficientIndex1"][0]).astype(int)
            idx2 = ak.to_numpy(chunk["events_EFTfitCoefficientIndex2"][0]).astype(int)
            wc_names = ["SM"] + decode_WCnames(ak.to_list(chunk["events_WCnames"][0]))
            meta = (idx1, idx2, wc_names)
            k00 = int(np.flatnonzero((idx1 == 0) & (idx2 == 0))[0])

        col = lambda name: ak.to_numpy(chunk[name])
        coefs = col("events_EFTfitCoefficients").astype(np.float64)
        coefs = np.column_stack([coefs, np.ones(len(coefs)), coefs[:, k00] ** 2])

        zh_pt, zh_m = col("events_ZH_pt"), col("events_ZH_M")
        sel = ((col("events_n_ak4jets") >= 5) & (col("events_n_b_outZH") == 2) &
               (col("events_ZH_bbvLscore") >= BBVL_WP) & (col("events_MET_pt") >= 20) &
               (zh_m >= 50) & (zh_m <= 200))
        genpt = col("events_genZHpt") if cfg["gen_split"] else None
        gmasks = gen_masks(label, cfg["gen_split"], genpt, len(sel))

        for region, (rmask, nn) in region_scores(col, model, cut).items():
            if region not in edges:
                continue
            for ch, (lo, hi) in CHANNELS.items():
                ch_sel = sel & rmask & (nn >= 0.0) & (zh_pt > lo) & (zh_pt <= hi)
                for proc, gmask in gmasks.items():
                    key, m = (ch, proc), ch_sel & gmask
                    for year, e in edges[region].items():
                        if key not in acc[region][year]:
                            acc[region][year][key] = np.zeros((len(e[ch]) - 1, len(MASS_BINS) - 1, coefs.shape[1]))
                        if m.any():
                            fill(acc[region][year][key], nn[m], zh_m[m], coefs[m], e[ch], MASS_BINS)
    return acc, meta


def merge_results(results):
    """Sum binned coefficients across inputs per (channel, process)."""
    merged, meta = {}, {}
    for acc, m in results:
        if m is None:
            continue
        for key, w in acc.items():
            if key in merged:
                prev = meta[key]
                if prev[2] != m[2] or not np.array_equal(prev[0], m[0]) or not np.array_equal(prev[1], m[1]):
                    raise ValueError(f"Coefficient layout differs between inputs for {key}")
                merged[key] += w
            else:
                merged[key], meta[key] = w.copy(), m
    return merged, meta


def read_inputs(model, cut, edges, processes=None, n_workers=4):
    """Read all EFT samples once: {region: {year: (merged, meta)}}."""
    processes = processes or PROCESSES
    jobs = [(label, cfg, path, model, cut, edges) for label, cfg in processes.items() for path in cfg["paths"]]
    if n_workers > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
            results = list(ex.map(process_path, *zip(*jobs)))
    else:
        results = [process_path(*job) for job in jobs]
    return {r: {y: merge_results([(acc[r][y], m) for acc, m in results]) for y in e} for r, e in edges.items()}


# ---------------- Building and checking the Combine entries (as EFTScaling.ipynb) ----------------
def combine_parameters(wc_names, wcs=None, wc_range=WC_RANGE):
    lo, hi = wc_range
    wanted = wc_names[1:] if wcs is None else wcs
    params = ["SM[1]"] + [f"{wc}[0,{lo},{hi}]" for wc in wanted if wc in wc_names]
    return params, [wc for wc in wanted if wc not in wc_names]


def build_entry(unrolled, idx1, idx2, wc_names, params, channel, process, log=print):
    sel = np.array([0 if p.split("[")[0] == "SM" else wc_names.index(p.split("[")[0]) for p in params])
    nbins, ncoef, n_nano = unrolled.shape[0], len(params), len(wc_names)
    M_full = np.zeros((nbins, n_nano, n_nano))           # coefficient of c_a*c_b = 2*M_ab for a != b
    for k, (a, b) in enumerate(zip(idx1, idx2)):
        if a == b:
            M_full[:, a, a] = unrolled[:, k]
        else:
            M_full[:, a, b] = M_full[:, b, a] = 0.5 * unrolled[:, k]
    M = M_full[:, sel][:, :, sel]
    sm = M[:, 0, 0].copy()
    good = sm > 0
    if np.any(~good):
        log(f"  [{channel}/{process}] {np.sum(~good)} bins with SM yield <= 0 -> pure SM scaling")
    M[good] /= sm[good][:, None, None]
    M[~good] = 0.0
    M[~good, 0, 0] = 1.0
    rows, cols = np.tril_indices(ncoef)
    entry = {"channel": channel, "process": process, "parameters": params, "scaling": M[:, rows, cols].tolist()}
    return entry, sm, good


def eval_scaling(entry, point):
    """Per-bin EFT/SM ratio of an entry at a WC point, e.g. {'ctgre': 0.5}."""
    names = [p.split("[")[0] for p in entry["parameters"]]
    theta = np.array([1.0 if n == "SM" else point.get(n, 0.0) for n in names])
    rows, cols = np.tril_indices(len(theta))
    return np.asarray(entry["scaling"]) @ (np.where(rows == cols, 1.0, 2.0) * theta[rows] * theta[cols])


def closure_check(entry, unrolled, idx1, idx2, wc_names, sm, good, point):
    """Direct polynomial sum vs. SM * (theta^T M theta); should agree to ~1e-12."""
    kept = {p.split("[")[0] for p in entry["parameters"]}
    c = np.array([1.0] + [point.get(n, 0.0) if n in kept else 0.0 for n in wc_names[1:]])
    direct = unrolled @ (c[idx1] * c[idx2])
    if not np.any(good):
        return 0.0
    predicted = sm * eval_scaling(entry, point)
    return float(np.max(np.abs(predicted[good] - direct[good]) / np.maximum(np.abs(direct[good]), 1e-30)))


def build_scaling(merged, meta, nn_edges, channel_name, min_events=20, max_rel_err=0.3, wcs=None,
                  wc_range=WC_RANGE, keep_diagonal=True, log=print):
    """interferenceModel entries of one region (fallback 'auto' of EFTScaling.ipynb).
    nn_edges: {channel: edges}; channel_name(ch) -> channel name in the combined card."""
    proc_totals = {}
    for (ch, proc), full in merged.items():
        ncoef = len(meta[(ch, proc)][0])
        proc_totals[proc] = proc_totals.get(proc, 0) + full[..., :ncoef].sum(axis=(0, 1))

    entries, worst = [], 0.0
    for (ch, proc) in sorted(merged):
        full = merged[(ch, proc)]
        idx1, idx2, wc_names = meta[(ch, proc)]
        ncoef = len(idx1)
        k00 = int(np.flatnonzero((idx1 == 0) & (idx2 == 0))[0])
        name = channel_name(ch)
        if not np.any(full[..., :ncoef]):
            log(f"  [{name}/{proc}] empty, skipped")
            continue
        params, missing = combine_parameters(wc_names, wcs, wc_range)
        if missing:
            log(f"  [{name}/{proc}] WCs not in sample (omitted): {missing}")

        unrolled_full = unroll(full, ch)
        unrolled = unrolled_full[:, :ncoef].copy()
        n_raw, sm2 = unrolled_full[:, ncoef].copy(), unrolled_full[:, ncoef + 1].copy()
        low_stat = unreliable(n_raw, unrolled[:, k00], sm2, min_events, max_rel_err)
        if keep_diagonal:
            low_stat &= ~diagonal_mask(ch, proc, nn_edges[ch])
        if np.any(low_stat):
            ch_tot = unrolled.sum(axis=0)
            ch_ok = not unreliable(np.array([n_raw.sum()]), np.array([ch_tot[k00]]),
                                   np.array([sm2.sum()]), min_events, max_rel_err)[0]
            unrolled[low_stat] = ch_tot if ch_ok else proc_totals[proc]   # normalization cancels in the ratio
            log(f"  [{name}/{proc}] {np.sum(low_stat)}/{len(low_stat)} bins smoothed -> "
                f"{'channel' if ch_ok else 'process'}-level ratio")

        entry, sm, good = build_entry(unrolled, idx1, idx2, wc_names, params, name, proc, log)
        rel = closure_check(entry, unrolled, idx1, idx2, wc_names, sm, good,
                            {p.split("[")[0]: 1.0 for p in params[1:]})
        worst = max(worst, rel)
        log(f"  [{name}/{proc}] bins={unrolled.shape[0]}, raw events={int(n_raw.sum())}, "
            f"SM sum={unrolled_full[:, k00].sum():.4g}, closure {rel:.1e} {'OK' if rel < 1e-8 else 'CHECK!'}")
        entries.append(entry)
    return entries, worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    scheme = ap.add_mutually_exclusive_group(required=True)
    scheme.add_argument("--cut", type=float, help="tt+LF/C score cut (makeDatacard.py --ttlf-cut)")
    scheme.add_argument("--argmax", action="store_true", help="argmax categories (makeDatacard.py --argmax)")
    ap.add_argument("--model", default="refactor_noOverlap_MultiWP", help="SPANet model of the scores")
    ap.add_argument("--edges", default=None,
                    help="--cut: SR edges JSON (default nn_edges_sigvsttbb_<model>_ttlf<cut>.json from makeDatacard.py)")
    ap.add_argument("--no-cr", action="store_true",
                    help="no CR / LF entries (cards made without --ttlf-cr / --argmax-lf)")
    ap.add_argument("--channel-format", default="{region}_{channel}",
                    help="channel names in the combined card (default SR_Zhpt1 ...)")
    ap.add_argument("--workers", type=int, default=4, help="1 = read serially")
    ap.add_argument("--all-wcs", action="store_true", help="every WC in the samples instead of DEFAULT_WCS")
    ap.add_argument("--min-events", type=int, default=20)
    ap.add_argument("--max-rel-err", type=float, default=0.3)
    args = ap.parse_args()

    if args.argmax:
        suffix, cut = "argmax", "argmax"
        edges = {"ZH": json.load(open(f"nn_edges_disc_{args.model}_argmaxZH.json")),
                 "TTB": json.load(open(f"nn_edges_ttbb_score_{args.model}_argmaxTTB.json"))}
        cr_name = "LF"
    else:
        suffix, cut = f"ttlf{args.cut:g}".replace(".", "p"), args.cut
        edges = {"SR": json.load(open(args.edges or f"nn_edges_sigvsttbb_{args.model}_{suffix}.json"))}
        cr_name = "CR"
    years = next(iter(edges.values()))
    if not args.no_cr:
        edges[cr_name] = {y: {ch: [0.0, 1.0] for ch in CHANNELS} for y in years}

    by_region = read_inputs(args.model, cut, edges, n_workers=args.workers)
    for year in sorted(years):
        entries, worst = [], 0.0
        for region in edges:
            print(f"--- {year} {region}")
            merged, meta = by_region[region][year]
            e, w = build_scaling(merged, meta, edges[region][year],
                                 lambda ch, r=region: args.channel_format.format(region=r, channel=ch),
                                 min_events=args.min_events, max_rel_err=args.max_rel_err,
                                 wcs=None if args.all_wcs else DEFAULT_WCS)
            entries += e
            worst = max(worst, w)
        out = f"combine_scaling_{suffix}_{year}.json"
        with open(out, "w") as f:
            json.dump(entries, f, indent=1)
        print(f"Wrote {len(entries)} entries ({', '.join(edges)}) to {out}, worst closure {worst:.1e}")


if __name__ == "__main__":
    main()
