#!/usr/bin/env python3
"""
Fit-free optimization of the Xbb working point, the discriminator and the NN binning,
directly from the analysis pickles.

Figure of merit (statistics + MC statistics), summed over all (NN, mass) bins of a channel:

    I_X = sum_bins S_X^2 / (T + V)      X = ttZ or ttH (measured signal, gen pT >= 200 GeV)
        T = total expected yield in the bin (everything else, incl. the other signal, is background)
        V = sum of w^2 of all MC in the bin  (penalizes bins with poor MC statistics)
    FOM = I_ttZ + I_ttH,   sigma(mu_X) ~ 1 / sqrt(I_X)

Tools (all print a table; validate the final choice with the full Combine fit):
    binning : greedy merging of a fine NN binning per channel, with an MC-statistics constraint
    disc    : tune disc = s / (s + a*ttbb + b*ttlf) per channel (a = b = 1 is the current disc)
    xbb     : compare Xbb thresholds (the pickles must contain events below the loosest one)

Usage:
    python3 optimize_binning.py binning --year 2018 --model overlapRemoval_MultiWP --out nn_edges_opt.json
    python3 optimize_binning.py disc    --year 2018 --model overlapRemoval_MultiWP
    python3 optimize_binning.py xbb     --year 2018 --model overlapRemoval_MultiWP --thresholds 0.80 0.85 0.9172 0.9734

Caveat: if the network was trained on the same MC events, optimize on the held-out events only
(otherwise the binning is tuned on overtrained scores).
"""
import argparse
import json

import numpy as np
import pandas as pd

PT_EDGES = [200, 300, 450, np.inf]
CHANNELS = ['Zhpt1', 'Zhpt2', 'Zhpt3']
MASS_EDGES = {'Zhpt1': [50, 80, 105, 200], 'Zhpt2': [50, 80, 105, 145, 200],
              'Zhpt3': [50, 80, 105, 145, 200]}
SIGNALS = ['ttZ', 'ttH']
BACKGROUNDS = ['TTBar', 'tt_B', 'SingleTop', 'ttX', 'VJets']


# ----------------------------------------------------------------------------- data
def select(df, year, xbb=None):
    """Analysis selection; xbb overrides the per-era Xbb working point."""
    from analysis_years import bbvl_wp
    if xbb is not None:
        thr = xbb
    elif 'era' in df:
        thr = df['era'].map({e: bbvl_wp(e) for e in df['era'].unique()})
    else:
        thr = bbvl_wp(year)
    m = ((df['n_ak4jets'] >= 5) & (df['n_b_outZH'] == 2) & (df['ZH_bbvLscore'] >= thr) &
         (df['MET_pt'] >= 20) & (df['ZH_M'] >= 50) & (df['ZH_M'] <= 200))
    return df[m]


def load_events(year, model, xbb=None, scale=1.0, pickle_dir=None):
    """One row per selected MC event: channel, mass, raw scores, weight, label."""
    from analysis_years import PICKLE_DIR
    from datacard_shapes import getZhbbWeight
    pickle_dir = pickle_dir or PICKLE_DIR
    frames = []
    for proc in SIGNALS + BACKGROUNDS:
        path = f"{pickle_dir.format(year=year)}/{proc}.pkl"
        try:
            df = pd.read_pickle(path)
        except FileNotFoundError:
            print(f"  [missing] {path}")
            continue
        df = select(df, year, xbb)
        out = pd.DataFrame({'pt': df['ZH_pt'].to_numpy(), 'M': df['ZH_M'].to_numpy(),
                            's': df[f'sig_score_{model}'].to_numpy(),
                            'ttbb': df[f'ttbb_score_{model}'].to_numpy(),
                            'ttlf': df[f'ttlf_score_{model}'].to_numpy(),
                            'w': scale * np.asarray(getZhbbWeight(df, year), dtype=float)})
        if proc in SIGNALS:
            out['label'] = np.where(df['genZHpt'].to_numpy() >= 200, proc, 'bkg')
        else:
            out['label'] = 'bkg'
        frames.append(out)
    ev = pd.concat(frames, ignore_index=True)
    ev['ch'] = np.digitize(ev['pt'], PT_EDGES) - 1
    return ev[ev['ch'] >= 0].reset_index(drop=True)


def disc(ev, a=1.0, b=1.0):
    den = ev['s'] + a * ev['ttbb'] + b * ev['ttlf']
    return np.where(den > 0, ev['s'] / np.where(den > 0, den, 1.0), 0.0)


# ----------------------------------------------------------------------------- figure of merit
def fine_edges(x, w, n):
    """Edges with equal measured-signal weight per bin, from 0 to 1."""
    order = np.argsort(x)
    xs, ws = np.asarray(x)[order], np.clip(np.asarray(w)[order], 0, None)
    cw = np.cumsum(ws)
    if cw[-1] <= 0:
        return np.linspace(0, 1, n + 1)
    e = np.interp(np.linspace(0, 1, n + 1)[1:-1], cw / cw[-1], xs)
    return np.unique(np.concatenate([[0.0], e, [1.0]]))


def cell_grid(d, ev, nn_edges, m_edges):
    """Per (fine NN bin, mass bin): S_ttZ, S_ttH, total T, MC variance V."""
    h = lambda wts: np.histogram2d(d, ev['M'], bins=[nn_edges, m_edges], weights=wts)[0]
    w = ev['w'].to_numpy()
    return {'ttZ': h(np.where(ev['label'] == 'ttZ', w, 0.0)),
            'ttH': h(np.where(ev['label'] == 'ttH', w, 0.0)),
            'T': h(w), 'V': h(w ** 2)}


def grouped(cells, bounds):
    """Sum fine NN rows into groups delimited by bounds (indices into the fine edges)."""
    out = {}
    for k, a in cells.items():
        c = np.vstack([np.zeros((1, a.shape[1])), np.cumsum(a, axis=0)])
        out[k] = c[bounds[1:]] - c[bounds[:-1]]
    return out


def information(g):
    T, V = g['T'], g['V']
    den = np.where(T + V > 0, T + V, 1.0)
    return {X: float(np.sum(np.where(T > 0, g[X] ** 2 / den, 0.0))) for X in SIGNALS}


def violations(g, max_rel):
    """Per NN group: True if any of its mass cells has no background with signal present,
    or a relative MC uncertainty of the total prediction above max_rel."""
    T, V = g['T'], g['V']
    S = g['ttZ'] + g['ttH']
    B = T - S
    bad = ((S > 0) & (B <= 0)) | ((T > 0) & (np.sqrt(V) > max_rel * T)) | ((T <= 0) & (S > 0))
    return bad.any(axis=1)


def optimize_edges(cells, edges, max_rel=0.3, max_bins=8):
    """Greedy merge of fine NN bins: first fix MC-statistics violations, then merge the
    pair costing the least information until at most max_bins NN bins remain."""
    bounds = np.arange(len(edges))
    fom = lambda b: sum(information(grouped(cells, b)).values())
    while len(bounds) > 2:
        bad = violations(grouped(cells, bounds), max_rel)
        if bad.any():
            i = int(np.flatnonzero(bad)[0])
            cands = [j for j in (i, i + 1) if 0 < j < len(bounds) - 1]      # internal bounds to drop
        elif len(bounds) - 1 > max_bins:
            cands = list(range(1, len(bounds) - 1))
        else:
            break
        best = max(cands, key=lambda j: fom(np.delete(bounds, j)))
        bounds = np.delete(bounds, best)
    g = grouped(cells, bounds)
    return np.asarray(edges)[bounds], information(g), bool(violations(g, max_rel).any())


def evaluate(ev, a=1.0, b=1.0, n_fine=40, max_rel=0.3, max_bins=8, verbose=False):
    """Optimized binning and information per channel for a given disc definition.
    a, b: scalars, or dicts per channel name."""
    res, tot = {}, {X: 0.0 for X in SIGNALS}
    for i, ch in enumerate(CHANNELS):
        sub = ev[ev['ch'] == i]
        if sub.empty:
            continue
        ai = a[ch] if isinstance(a, dict) else a
        bi = b[ch] if isinstance(b, dict) else b
        d = disc(sub, ai, bi)
        meas = sub['label'].isin(SIGNALS).to_numpy()
        fe = fine_edges(d[meas], sub['w'].to_numpy()[meas], n_fine)
        cells = cell_grid(d, sub, fe, MASS_EDGES[ch])
        edges, info, still_bad = optimize_edges(cells, fe, max_rel, max_bins)
        res[ch] = {'edges': [float(e) for e in edges], 'info': info, 'mc_ok': not still_bad}
        for X in SIGNALS:
            tot[X] += info[X]
        if verbose:
            print(f"  {ch}: {len(edges) - 1} NN bins, I_ttZ={info['ttZ']:.3f}, I_ttH={info['ttH']:.3f}"
                  f"{'' if not still_bad else '  (MC constraint not met)'}\n      edges {np.round(edges, 4).tolist()}")
    return res, tot


def evaluate_equal_bins(ev, a=1.0, b=1.0, n_fine=40, max_rel=0.3, max_bins=8, verbose=False):
    """Like evaluate, but every channel ends with the same number of NN bins
    (make_datacard.py stacks the channels into one array, so they must match)."""
    n = max_bins
    while True:
        res, tot = evaluate(ev, a, b, n_fine, max_rel, n)
        counts = {len(r['edges']) - 1 for r in res.values()}
        if len(counts) <= 1 or n <= 1:
            break
        n = min(counts)
    if verbose:
        evaluate(ev, a, b, n_fine, max_rel, n, verbose=True)
    return res, tot


def sigma_line(tot):
    return "  ".join(f"sigma({X})~{1 / np.sqrt(v):.3f}" if v > 0 else f"sigma({X})=inf"
                     for X, v in tot.items())


# ----------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['binning', 'disc', 'xbb'])
    ap.add_argument('--year', default='2018')
    ap.add_argument('--model', default='overlapRemoval_MultiWP')
    ap.add_argument('--scale', type=float, default=2.31, help='luminosity scale (2018 -> 138 fb^-1: 2.31)')
    ap.add_argument('--max-bins', type=int, default=8, help='max NN bins per channel')
    ap.add_argument('--max-rel', type=float, default=0.3, help='max relative MC uncertainty per bin')
    ap.add_argument('--n-fine', type=int, default=40)
    ap.add_argument('--out', default=None, help='binning mode: write optimized edges (JSON)')
    ap.add_argument('--thresholds', type=float, nargs='+', default=[0.85, 0.9172, 0.9734, 0.9880])
    ap.add_argument('--grid', type=float, nargs='+', default=[0.25, 0.5, 1.0, 2.0, 4.0])
    args = ap.parse_args()
    kw = dict(n_fine=args.n_fine, max_rel=args.max_rel, max_bins=args.max_bins)

    if args.mode == 'binning':
        ev = load_events(args.year, args.model, scale=args.scale)
        res, tot = evaluate_equal_bins(ev, verbose=True, **kw)
        print(f"\nTotal: I_ttZ={tot['ttZ']:.3f}  I_ttH={tot['ttH']:.3f}   {sigma_line(tot)}")
        if args.out:
            with open(args.out, 'w') as f:
                json.dump({args.year: {ch: r['edges'] for ch, r in res.items()}}, f, indent=1)
            print(f"Wrote {args.out}  (use as DataCardShapes(..., edges_file=...) and as NN_BINS)")

    elif args.mode == 'disc':
        ev = load_events(args.year, args.model, scale=args.scale)
        _, ref = evaluate(ev, **kw)
        print(f"reference a=b=1: FOM={sum(ref.values()):.3f}   {sigma_line(ref)}")
        best = {}
        for i, ch in enumerate(CHANNELS):
            sub = ev[ev['ch'] == i]
            scores = {(a, b): sum(evaluate(sub.assign(ch=i), a, b, **kw)[1].values())
                      for a in args.grid for b in args.grid}
            (ba, bb), fom = max(scores.items(), key=lambda kv: kv[1])
            best[ch] = (ba, bb)
            print(f"  {ch}: best a(ttbb)={ba}, b(ttlf)={bb}: FOM {scores[(1.0, 1.0)]:.3f} -> {fom:.3f}")
        a = {ch: v[0] for ch, v in best.items()}
        b = {ch: v[1] for ch, v in best.items()}
        _, tot = evaluate(ev, a, b, **kw)
        print(f"with per-channel best weights: FOM={sum(tot.values()):.3f}   {sigma_line(tot)}")

    elif args.mode == 'xbb':
        print(f"{'Xbb cut':>8} {'I_ttZ':>8} {'I_ttH':>8} {'FOM':>8}   estimated precision")
        for thr in args.thresholds:
            ev = load_events(args.year, args.model, xbb=thr, scale=args.scale)
            _, tot = evaluate(ev, **kw)
            print(f"{thr:8.4f} {tot['ttZ']:8.3f} {tot['ttH']:8.3f} {sum(tot.values()):8.3f}   {sigma_line(tot)}")


if __name__ == '__main__':
    main()