#!/usr/bin/env python3
"""
b-tag shape scale-factor normalization (BTV recommendation).

The shape-method b-tag weights must not change the expected yield before any b-tag
requirement. For every era, sample (process / flavour split, e.g. 'TTToSemiLeptonic__tt+LF')
and AK4 jet multiplicity bin, and separately for the nominal SF and each SF variation:

    r = sum(w_before) / sum(w_after)
      w_before = event weight without the b-tag SF
      w_after  = w_before * b-tag SF (nominal or varied)

computed in a selection WITHOUT any AK4 b-tag requirement. r is then multiplied into the
b-tag SF columns of the analysis pickles (coffea_to_pkl.py --btag-norm btag_norm.json).

Inputs: a .coffea output with an extra category (default 'btag_calib') that has the analysis
preselection minus all AK4 b-tag requirements, exporting the columns
    n_ak4jets, genWeight, the SF columns in WEIGHT_COLUMNS, btag_sf and btag_sf<src>_up/_down

Usage:
    python3 btag_norm.py -i output_all.coffea -o btag_norm.json [--category btag_calib]
"""
import argparse
import json

import numpy as np
import pandas as pd

# Jet multiplicity bins: [lo, hi). The last bin collects everything above.
NJET_EDGES = [0, 5, 6, 7, 8, 9, 1000]
BTAG_SOURCES = ['hf', 'lf', 'cferr1', 'cferr2', 'hfstats1', 'hfstats2', 'lfstats1', 'lfstats2']
VARIATIONS = {'nominal': 'btag_sf',
              **{f'{s}_{d}': f'btag_sf{s}_{d}' for s in BTAG_SOURCES for d in ('up', 'down')}}
MIN_EVENTS = 200          # raw MC events per jet bin, else fall back to the sample-inclusive r

# Weight factors other than the b-tag SF (same list as datacard_shapes.WEIGHT_COLUMNS)
OTHER_WEIGHTS = ['topptWeight', 'ele_reco_sf', 'ele_id_sf', 'ele_trig_sf',
                 'mu_id_sf', 'mu_iso_sf', 'mu_reco_sf', 'mu_trig_sf', 'xbb_sf', 'puWeight']


def njet_bin(njets):
    return np.clip(np.digitize(np.asarray(njets), NJET_EDGES) - 1, 0, len(NJET_EDGES) - 2)


def weight_without_btag(df):
    w = df['norm_weight'].to_numpy() * np.sign(df['genWeight'].to_numpy())
    for col in OTHER_WEIGHTS:
        if col in df:
            w = w * df[col].to_numpy()
    return w


def compute_ratios(df):
    """r per variation and jet bin for ONE (era, sample). Returns {variation: {...}}."""
    w0 = weight_without_btag(df)
    idx = njet_bin(df['n_ak4jets'])
    nb = len(NJET_EDGES) - 1
    n_raw = np.bincount(idx, minlength=nb)
    before = np.bincount(idx, weights=w0, minlength=nb)
    out = {}
    for var, col in VARIATIONS.items():
        if col not in df:
            continue
        sf = df[col].to_numpy()
        after = np.bincount(idx, weights=w0 * sf, minlength=nb)
        r_incl = before.sum() / after.sum() if after.sum() > 0 else 1.0
        ok = (n_raw >= MIN_EVENTS) & (after > 0) & (before > 0)
        r = np.where(ok, before / np.where(after > 0, after, 1.0), r_incl)
        out[var] = {'r': [float(v) for v in r], 'r_inclusive': float(r_incl),
                    'fallback_bins': [int(i) for i in np.flatnonzero(~ok)]}
    out['_n_raw'] = [int(n) for n in n_raw]
    return out


def read_calib_category(coffea_file, category, samples=None):
    """{(era, sample): DataFrame} from the calibration category of a .coffea output."""
    from coffea.util import load
    from analysis_years import find_era
    from outputToPkl import get_norm

    filein = load(coffea_file)
    genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')
    needed = ['n_ak4jets', 'genWeight'] + OTHER_WEIGHTS + list(VARIATIONS.values())
    frames = {}
    for sample, datasets in filein['columns'].items():
        if samples and sample not in samples:
            continue
        for dataset, cats in datasets.items():
            if category not in cats:
                continue
            cols = cats[category]['nominal']
            era = find_era(dataset)
            data = {c: cols[f'events_{c}'].value for c in needed if f'events_{c}' in cols}
            if 'n_ak4jets' not in data or 'btag_sf' not in data:
                print(f"  [skip] {sample}/{dataset}: n_ak4jets or btag_sf missing in '{category}'")
                continue
            df = pd.DataFrame(data)
            df['norm_weight'] = get_norm(sample, dataset, genweight_df, era)
            key = (era, sample)
            frames[key] = pd.concat([frames[key], df], ignore_index=True) if key in frames else df
    return frames


def build_table(frames):
    table = {'njet_edges': NJET_EDGES}
    for (era, sample), df in sorted(frames.items()):
        res = compute_ratios(df)
        table.setdefault(era, {})[sample] = res
        nom = res.get('nominal', {})
        print(f"  {era:13} {sample:32} r(nominal) per jet bin: "
              + " ".join(f"{v:.3f}" for v in nom.get('r', []))
              + (f"   fallback bins {nom['fallback_bins']}" if nom.get('fallback_bins') else ""))
    return table


def apply_btag_norm(df, era, sample, table):
    """Multiply the b-tag SF columns of an analysis DataFrame by r(era, sample, njets).
    Variations get their own r, so btag_sf<src>_up / btag_sf keeps the BTV-normalized ratio."""
    entry = table.get(era, {}).get(sample)
    if entry is None:
        print(f"  [WARNING] no b-tag normalization for {era}/{sample}; left unnormalized")
        return df
    idx = njet_bin(df['n_ak4jets'])
    for var, col in VARIATIONS.items():
        if col in df and var in entry:
            df[col] = df[col].to_numpy() * np.asarray(entry[var]['r'])[idx]
    return df


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', '-i', required=True, help='.coffea output containing the calibration category')
    ap.add_argument('--output', '-o', default='btag_norm.json')
    ap.add_argument('--category', default='btag_calib')
    args = ap.parse_args()
    table = build_table(read_calib_category(args.input, args.category))
    with open(args.output, 'w') as f:
        json.dump(table, f, indent=1)
    print(f"Wrote {args.output}")