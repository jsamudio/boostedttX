#!/usr/bin/env python3
"""
Xbb working-point scan of the Run 3 analysis, read directly from the analysis output (no pickles yet).

Events: btag_mask category of the main inference (signal, bkg and minor groups), with the selection of
the pickles (as in Run 2: n_ak4jets >= 5, n_b_outZH == 2, MET >= 20 GeV, 50 <= ZH_M <= 200 GeV) and a
cut on the ZH candidate score ZH_bbvLscore (GloParT Xbb / (Xbb + QCD) of the leading-Xbb fat jet).

Figure of merit (as optimize_binning.py --fom fisher in Run 2): the measured processes (--pois, default
ttZ ttH tt_B) are fitted together, with the Fisher matrix over all bins

    F_ab = sum_bins S_a S_b / (T + V),   I_X = 1 / (F^-1)_XX,   sigma(mu_X) ~ 1 / sqrt(I_X)

    T = total expected yield in the bin, V = sum of w^2 (MC statistics). Bins: pT channel (200, 300,
    450, inf) x ZH mass ([50, 80, 105, 200] in Zhpt1, [50, 80, 105, 145, 200] above) [x Xbb category].
    ttZ / ttH: signal with gen Z/H pT >= 200 GeV (the rest is background); tt_B: all tt+B events.

Modes:
    cut     : one threshold per value of --cuts. Prints the efficiencies, the yields and sigma(mu) of each
              POI, also relative to the loosest cut.
    multiwp : preselection at --presel, events split into Xbb categories [t0, t1), [t1, t2), ..., [tn, 1].
              All combinations of --n-cat - 1 thresholds from --grid (default: the 10%..90% quantiles of
              the signal and of the rest, above the preselection) are ranked by FOM = mean_X I_X /
              I_X(continuous), where "continuous" bins the score in --n-fine bins (half signal, half
              background quantiles): close to the upper bound of what the score can add. This is the
              proxy for the SPANet Xbb input as WP categories instead of the continuous score.

    composition : process composition after the selection vs the cut: yields of ttZ / ttH (gen pT >=
              200 GeV), sig_lowpt, tt_B, tt+C, tt+LF, minor, total and per pT channel, and the ratios
              tt+LF / S, tt+C / S, tt_B / S (S = ttZ + ttH, or ttH with --signal ttH). With --ref
              <Run2/xbb_reference_<year>.json> (Run 2 after Xbb >= LP): the Run 3 cut at which each Run 3
              ratio equals the Run 2 one, with the efficiencies of each class there (relative to no cut).
              The cut matching tt+LF / S is the Run 3 cut with the Run 2 light rejection relative to the
              signal (b tagger, sqrt(s) and Xbb SFs also differ: see Run2/xbb_reference.py).

Usage (inside the image):
    python3 scan_xbb_wp.py cut -i outputs/nominal/{signal,bkg,minor}/merged/*.coffea   (PBS per-dataset merge; else output_all.coffea)
    python3 scan_xbb_wp.py cut -i ... --cuts 0.1 0.3 0.5 0.7 0.9 --pois ttZ ttH
    python3 scan_xbb_wp.py multiwp -i ... --presel 0.3 --n-cat 3
    python3 scan_xbb_wp.py composition -i ... --ref ../Run2/xbb_reference_2018.json [--signal ttH]

Not included: systematics other than tt_B floating (the Xbb SF uncertainty grows with a tighter WP and
is not measured yet), the SPANet (the optimum moves once the network exists, since it uses the Xbb
information; the multiwp mode only uses the candidate fat jet, while the network sees every fat jet),
and process normalizations beyond xsec x lumi.
"""
import argparse
import itertools

import numpy as np

import samples
from outvars import nominal_weight_factors

CATEGORY = "btag_mask"
VARIATION = "nominal"
PT_EDGES = [200, 300, 450, np.inf]
MASS_EDGES = [[50, 80, 105, 200], [50, 80, 105, 145, 200], [50, 80, 105, 145, 200]]
SIG_GENPT = 200.0
DEFAULT_CUTS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98]
# composition mode: also the light-mistag cuts of derive_xbb_wps.py (10/5/1%) and the multiwp best set
COMP_CUTS = [0.0, 0.0026, 0.0068, 0.0542, 0.1, 0.2, 0.3, 0.5, 0.7, 0.877, 0.9, 0.95, 0.975]

# Output keys -> process (as in the Run 2 pickles: tt+B from the TTbb samples, tt+LF / tt+C from the
# inclusive ttbar samples)
SIGNAL_KEYS = {f"{s}__{m}": samples.PROCESS_LABELS[s]
               for s in samples.SIGNAL_SAMPLES for m in samples.GENMATCH_SPLIT}
BKG_KEYS = {f"{s}__tt+B": "tt_B" for s in samples.TTBAR_SAMPLES if s.startswith("TTbb")}
BKG_KEYS.update({f"{s}__{c}": "TTBar" for s in samples.TTBAR_SAMPLES if s.startswith("TTTo")
                 for c in ("tt+LF", "tt+C")})
BKG_KEYS.update({s: samples.PROCESS_LABELS[s] for s in samples.MINOR_SAMPLES + samples.QCD_SAMPLES})

READ = ["ZH_pt", "ZH_M", "ZH_bbvLscore", "n_ak4jets", "n_b_outZH", "MET_pt", "genWeight", "norm_weight"]


# ----------------------------------------------------------------------------- data
def load_events(paths, pois):
    """One array per field over all selected MC events (any Xbb score): pt, M, xbb, w, proc, poi, bin."""
    from coffea.util import load
    out = {k: [] for k in ("pt", "M", "xbb", "w", "proc", "poi", "cls")}
    missing_factors = set()
    for path in paths:
        print(f"Reading {path}")
        filein = load(path)
        gw = filein["sum_signOf_genweights"]
        for key, proc in {**SIGNAL_KEYS, **BKG_KEYS}.items():
            for dataset, cats in filein["columns"].get(key, {}).items():
                cols = cats.get(CATEGORY, {}).get(VARIATION)
                if cols is None:
                    continue
                c = {v: np.asarray(cols[f"events_{v}"].value) for v in READ}
                w = c["norm_weight"] * np.sign(c["genWeight"]) / gw[dataset]
                for f in nominal_weight_factors:
                    if f"events_{f}" in cols:
                        w = w * np.asarray(cols[f"events_{f}"].value)
                    else:
                        missing_factors.add(f)
                sel = ((c["n_ak4jets"] >= 5) & (c["n_b_outZH"] == 2) & (c["MET_pt"] >= 20)
                       & (c["ZH_M"] >= 50) & (c["ZH_M"] <= 200) & (c["ZH_pt"] >= PT_EDGES[0]))
                if key in SIGNAL_KEYS:
                    genpt = np.asarray(cols["events_genZHpt"].value)
                    poi = np.where(genpt >= SIG_GENPT, proc, "bkg") if proc in pois else np.full(len(w), "bkg")
                    cls = np.where(genpt >= SIG_GENPT, proc, "sig_lowpt")
                else:
                    poi = np.full(len(w), proc if proc in pois else "bkg")
                    cls = np.full(len(w), composition_class(key, proc))
                out["pt"].append(c["ZH_pt"][sel])
                out["M"].append(c["ZH_M"][sel])
                out["xbb"].append(c["ZH_bbvLscore"][sel])
                out["w"].append(w[sel])
                out["proc"].append(np.full(int(sel.sum()), proc))
                out["poi"].append(poi[sel])
                out["cls"].append(cls[sel])
        del filein
    if not out["w"]:
        raise RuntimeError(f"No events of the {CATEGORY} category in {paths}")
    if missing_factors:
        print(f"  [warning] weight factors not in the output (taken as 1): {sorted(missing_factors)}")
    ev = {k: np.concatenate(v) for k, v in out.items()}
    ev["bin"] = mass_pt_bin(ev["pt"], ev["M"])
    return ev


def composition_class(key, proc):
    """Background class of an output key, as Run2/xbb_reference.py: tt_B, tt+C, tt+LF or minor."""
    if proc == "tt_B":
        return "tt_B"
    if proc == "TTBar":
        return "tt+C" if key.endswith("__tt+C") else "tt+LF"
    return "minor"


def mass_pt_bin(pt, M):
    """Index of the (pT channel, mass bin) cell, consecutive over the channels."""
    ch = np.digitize(pt, PT_EDGES) - 1
    out = np.full(len(pt), -1)
    offset = 0
    for i, edges in enumerate(MASS_EDGES):
        m = ch == i
        out[m] = offset + np.clip(np.digitize(M[m], edges) - 1, 0, len(edges) - 2)
        offset += len(edges) - 1
    return out


N_MASS_PT = sum(len(e) - 1 for e in MASS_EDGES)


def subset(ev, mask):
    return {k: v[mask] for k, v in ev.items()}


# ----------------------------------------------------------------------------- figure of merit
def information(ev, pois, cat=None, ncat=1):
    """I_X = 1 / (F^-1)_XX per POI, over the (pT, mass[, Xbb category]) bins of the events in ev."""
    idx = ev["bin"] if cat is None else ev["bin"] * ncat + cat
    n = N_MASS_PT * ncat
    w = ev["w"]
    S = {X: np.bincount(idx, np.where(ev["poi"] == X, w, 0.0), n) for X in pois}
    T = np.bincount(idx, w, n)
    V = np.bincount(idx, w ** 2, n)
    ok = T > 0
    den = np.where(ok, T + V, 1.0)
    F = np.array([[float(np.sum(np.where(ok, S[a] * S[b] / den, 0.0))) for b in pois] for a in pois])
    info = {}
    for i, X in enumerate(pois):
        r = [j for j in range(len(pois)) if j != i]
        v = F[i, i] - (F[i, r] @ np.linalg.pinv(F[np.ix_(r, r)]) @ F[r, i] if r else 0.0)
        info[X] = float(max(v, 0.0))
    return info


def sigma(info):
    return {X: (1 / np.sqrt(I) if I > 0 else np.inf) for X, I in info.items()}


def yields(ev, pois):
    out = {X: float(ev["w"][ev["poi"] == X].sum()) for X in pois}
    out["total"] = float(ev["w"].sum())
    out["TTBar"] = float(ev["w"][ev["proc"] == "TTBar"].sum())
    return out


def quantile_edges(x, w, qs):
    """Weighted quantiles of x (weights clipped at 0)."""
    order = np.argsort(x)
    xs, ws = x[order], np.clip(w[order], 0, None)
    cw = np.cumsum(ws)
    return np.interp(qs, cw / cw[-1], xs)


def balanced_edges(ev, pois, n):
    """n - 1 quantiles of the measured signal (ttZ + ttH) and n - 1 of the rest (backgrounds and tt_B),
    merged: thresholds both where the signal falls off and where the backgrounds do."""
    sig = np.isin(ev["poi"], [X for X in pois if X in ("ttZ", "ttH")])
    qs = np.linspace(0, 1, n + 1)[1:-1]
    edges = [quantile_edges(ev["xbb"][m], ev["w"][m], qs) for m in (sig, ~sig) if m.any()]
    return np.unique(np.round(np.concatenate(edges), 3))


def xbb_category(x, thresholds):
    return np.digitize(x, thresholds)


# ----------------------------------------------------------------------------- modes
def run_cut(ev, pois, cuts):
    cuts = sorted(cuts)
    base = yields(subset(ev, ev["xbb"] >= cuts[0]), pois)
    ref = None
    head = (f"{'cut':>6} " + " ".join(f"{'eff ' + X:>10}" for X in pois) + f" {'eff TTBar':>10} "
            + " ".join(f"{'N ' + X:>9}" for X in pois) + f" {'N total':>10} "
            + " ".join(f"{'sig ' + X:>9}" for X in pois) + "  " + " ".join(f"{'rel ' + X:>8}" for X in pois))
    print(f"\nXbb cut scan on ZH_bbvLscore (efficiencies relative to cut {cuts[0]}); POIs fitted together: {pois}")
    print(head)
    for cut in cuts:
        sub = subset(ev, ev["xbb"] >= cut)
        y = yields(sub, pois)
        s = sigma(information(sub, pois))
        ref = ref or s
        eff = lambda k: y[k] / base[k] if base[k] > 0 else 0.0
        print(f"{cut:6.3f} " + " ".join(f"{eff(X):10.3f}" for X in pois) + f" {eff('TTBar'):10.3f} "
              + " ".join(f"{y[X]:9.1f}" for X in pois) + f" {y['total']:10.1f} "
              + " ".join(f"{s[X]:9.3f}" for X in pois) + "  "
              + " ".join(f"{s[X] / ref[X]:8.3f}" for X in pois))
    print("sig X = expected sigma(mu_X), stat + MC stat only; rel = sigma / sigma at the loosest cut")


COMP_CLASSES = ["ttZ", "ttH", "sig_lowpt", "tt_B", "tt+C", "tt+LF", "minor"]
RATIOS = {"tt+LF": "LF", "tt+C": "C", "tt_B": "B"}


def yields_above(ev, cuts):
    """{class: yields of the events with xbb >= cut, per cut}, total and per pT channel."""
    order = np.argsort(ev["xbb"])
    x = ev["xbb"][order]
    pos = np.searchsorted(x, cuts, side="left")
    ch = np.digitize(ev["pt"][order], PT_EDGES) - 1

    def above(mask):
        w = np.where(mask, ev["w"][order], 0.0)
        tail = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])     # tail[i] = sum of w[i:]
        return tail[pos]
    tot = {c: above(ev["cls"][order] == c) for c in COMP_CLASSES}
    chans = [{c: above((ev["cls"][order] == c) & (ch == i)) for c in COMP_CLASSES} for i in range(len(PT_EDGES) - 1)]
    return tot, chans


def ratios(y, signal):
    s = y["ttH"] if signal == "ttH" else y["ttZ"] + y["ttH"]
    return {c: np.divide(y[c], s, out=np.full_like(np.asarray(s, dtype=float), np.nan), where=np.asarray(s) > 0)
            for c in RATIOS}


def composition_table(rows, signal, title):
    print(f"\n{title}")
    sname = "ttH" if signal == "ttH" else "S"
    print(f"{'':>13} " + " ".join(f"{c:>10}" for c in COMP_CLASSES)
          + " " + " ".join(f"{v + '/' + sname:>8}" for v in RATIOS.values()))
    for name, y in rows:
        r = ratios(y, signal)
        print(f"{name:>13} " + " ".join(f"{float(y[c]):10.1f}" for c in COMP_CLASSES)
              + " " + " ".join(f"{float(r[c]):8.2f}" for c in RATIOS))


def run_composition(ev, cuts, ref_path, signal):
    cuts = sorted(cuts)
    tot, chans = yields_above(ev, np.array(cuts))
    at = lambda y, i: {c: y[c][i] for c in COMP_CLASSES}
    composition_table([(f"cut {c:.4f}", at(tot, i)) for i, c in enumerate(cuts)], signal,
                      "Run 3 composition after the selection vs the cut on ZH_bbvLscore (all pT channels)")
    for k in range(len(PT_EDGES) - 1):
        composition_table([(f"cut {c:.4f}", at(chans[k], i)) for i, c in enumerate(cuts)], signal,
                          f"pT channel {PT_EDGES[k]:.0f}-{PT_EDGES[k + 1]:.0f} GeV")
    if ref_path is None:
        return

    import json
    ref = json.load(open(ref_path))
    rows = [("total", ref["total"])] + [(f"{ref['pt_edges'][k]:.0f}-{ref['pt_edges'][k + 1] or np.inf:.0f}", c)
                                         for k, c in enumerate(ref["channels"])]
    composition_table(rows, signal, f"Run 2 reference {ref['year']} (after Xbb >= LP), {ref_path}")

    # fine grid: every 0.1% quantile of the score plus a linear grid
    grid = np.unique(np.concatenate([np.quantile(ev["xbb"], np.linspace(0, 1, 1001)), np.linspace(0, 1, 1001)]))
    gtot, gch = yields_above(ev, grid)
    none_ = {c: float(ev["w"][ev["cls"] == c].sum()) for c in COMP_CLASSES}      # no Xbb cut
    print(f"\nRun 3 cut matching each Run 2 ratio ({'ttH' if signal == 'ttH' else 'S = ttZ + ttH'}, all pT "
          f"channels); efficiencies relative to no Xbb cut:")
    print(f"{'ratio':>8} {'Run 2':>7} {'cut':>8} " + " ".join(f"{'eff ' + c:>11}" for c in COMP_CLASSES)
          + "   per-channel ratio Run 3 / Run 2")
    rref = ratios(ref["total"], signal)
    for c, short in RATIOS.items():
        r3 = ratios(gtot, signal)[c]
        ok = np.where(np.isfinite(r3) & (r3 <= rref[c]))[0]
        if not len(ok):
            print(f"{short + '/S':>8} {rref[c]:7.2f}   not reached (Run 3 min {np.nanmin(r3):.2f})")
            continue
        i = ok[0]
        eff = " ".join(f"{gtot[k][i] / none_[k] if none_[k] > 0 else np.nan:11.3f}" for k in COMP_CLASSES)
        per = []
        for j in range(len(PT_EDGES) - 1):
            r3j = float(ratios({k: gch[j][k][i] for k in COMP_CLASSES}, signal)[c])
            r2j = float(ratios(ref["channels"][j], signal)[c])
            per.append(f"{r3j / r2j:.2f}" if r2j > 0 else "n/a")
        per = " ".join(per)
        print(f"{short + '/S':>8} {rref[c]:7.2f} {grid[i]:8.4f} {eff}   {per}")
    print("(first cut from below where the Run 3 ratio is at or below the Run 2 one; the ratios need not be "
          "monotonic in the cut)")


def run_multiwp(ev, pois, presel, n_cat, grid, n_fine, top):
    sub = subset(ev, ev["xbb"] >= presel)
    if grid is None:
        grid = [g for g in balanced_edges(sub, pois, 10) if g > presel]
    grid = sorted(float(g) for g in grid)
    print(f"\nPreselection ZH_bbvLscore >= {presel}: {len(sub['w'])} MC events; POIs fitted together: {pois}")
    print(f"Threshold grid ({len(grid)}): {[round(g, 4) for g in grid]}")

    one = information(sub, pois)
    fine = balanced_edges(sub, pois, n_fine // 2)
    cont = information(sub, pois, xbb_category(sub["xbb"], fine), n_fine)
    gain = lambda info: float(np.mean([info[X] / cont[X] for X in pois if cont[X] > 0]))
    print(f"\n{'scheme':<34} " + " ".join(f"{'sig ' + X:>9}" for X in pois) + f" {'FOM':>7}")
    line = lambda name, info: print(f"{name:<34} " + " ".join(f"{s:9.3f}" for s in sigma(info).values())
                                    + f" {gain(info):7.3f}")
    line(f"continuous ({n_fine} quantile bins)", cont)
    line("single cut (1 category)", one)

    results = []
    for th in itertools.combinations(grid, n_cat - 1):
        info = information(sub, pois, xbb_category(sub["xbb"], list(th)), n_cat)
        results.append((gain(info), th, info))
    results.sort(key=lambda r: -r[0])
    print(f"\nBest {min(top, len(results))} of {len(results)} sets of {n_cat - 1} thresholds above {presel} "
          f"(FOM = mean_X I_X / I_X(continuous)):")
    for fom, th, info in results[:top]:
        line("[" + ", ".join(f"{t:.4f}" for t in th) + "]", info)

    best = results[0][1]
    cat = xbb_category(sub["xbb"], list(best))
    edges = [presel] + list(best) + [1.0]
    print(f"\nYields per category of the best set:")
    print(f"{'category':<20} " + " ".join(f"{X:>9}" for X in pois) + f" {'TTBar':>9} {'total':>10}")
    for k in range(n_cat):
        y = yields(subset(sub, cat == k), pois)
        name = f"[{edges[k]:.4f}, {edges[k + 1]:.4f}" + ("]" if k == n_cat - 1 else ")")
        print(f"{name:<20} " + " ".join(f"{y[X]:9.1f}" for X in pois) + f" {y['TTBar']:9.1f} {y['total']:10.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["cut", "multiwp", "composition"])
    ap.add_argument("-i", "--inputs", nargs="+", required=True, help=".coffea files of the signal, bkg, minor groups (merged/*.coffea or output_all.coffea)")
    ap.add_argument("--pois", nargs="+", default=["ttZ", "ttH", "tt_B"])
    ap.add_argument("--cuts", nargs="+", type=float, default=DEFAULT_CUTS, help="cut mode: thresholds")
    ap.add_argument("--presel", type=float, default=0.0, help="multiwp mode: preselection threshold")
    ap.add_argument("--n-cat", type=int, default=3, help="multiwp mode: number of Xbb categories above --presel")
    ap.add_argument("--grid", nargs="+", type=float, default=None,
                    help="multiwp mode: candidate thresholds (default: signal and background deciles above --presel)")
    ap.add_argument("--n-fine", type=int, default=20, help="multiwp mode: quantile bins of the continuous reference")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--ref", help="composition mode: Run 2 reference JSON (Run2/xbb_reference.py)")
    ap.add_argument("--signal", choices=["ttZ+ttH", "ttH"], default="ttZ+ttH",
                    help="composition mode: signal of the ratios (ttH: independent of the ttZ cross-section update)")
    args = ap.parse_args()

    ev = load_events(args.inputs, args.pois)
    ev = subset(ev, ev["bin"] >= 0)
    if args.mode == "cut":
        run_cut(ev, args.pois, args.cuts)
    elif args.mode == "composition":
        cuts = args.cuts if args.cuts is not DEFAULT_CUTS else COMP_CUTS
        run_composition(ev, cuts, args.ref, args.signal)
    else:
        if args.n_cat < 2:
            ap.error("--n-cat must be >= 2")
        run_multiwp(ev, args.pois, args.presel, args.n_cat, args.grid, args.n_fine, args.top)


if __name__ == "__main__":
    main()
