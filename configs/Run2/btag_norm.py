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
b-tag SF columns of the analysis pickles (outputToPkl.py --btag-norm btag_norm.json).

Inputs: one or more .coffea outputs of the inference configs. The processor
(spanet_inference.py) accumulates the per-jet-bin sums directly in output['btag_calib'],
so no per-event columns are needed. Older outputs that instead have a 'btag_calib'
column category are still supported (--category).

r is derived separately for every shape variation found in the inputs: for AK4 JES
variations btag_sf is the JES-varied SF, so it needs its own normalization. The JSON keeps
the nominal table at the top level and the others under 'shape_variations'.

Usage:
    python3 btag_norm.py -i outputs/nominal/*/output_all.coffea outputs/JER/*/output_all.coffea -o btag_norm.json
"""
import argparse
import json

import numpy as np
import pandas as pd

from outvars import nominal_weight_factors

# Jet multiplicity bins: [lo, hi). The last bin collects everything above.
NJET_EDGES = [0, 5, 6, 7, 8, 9, 1000]
BTAG_SOURCES = ['hf', 'lf', 'cferr1', 'cferr2', 'hfstats1', 'hfstats2', 'lfstats1', 'lfstats2']
VARIATIONS = {'nominal': 'btag_sf',
              **{f'{s}_{d}': f'btag_sf{s}_{d}' for s in BTAG_SOURCES for d in ('up', 'down')}}
MIN_EVENTS = 200          # raw MC events per jet bin, else fall back to the sample-inclusive r

# Weight factors other than the b-tag SF
OTHER_WEIGHTS = [w for w in nominal_weight_factors if w != 'btag_sf']


def njet_bin(njets):
    return np.clip(np.digitize(np.asarray(njets), NJET_EDGES) - 1, 0, len(NJET_EDGES) - 2)


def calib_sums(njets, w0, sfs):
    """Per-jet-bin sums for one set of events (used by the processor and for old outputs).
    sfs: {variation: per-event b-tag SF}. All sums are without the per-dataset normalization."""
    idx = njet_bin(njets)
    nb = len(NJET_EDGES) - 1

    def sums(weights=None):
        # Always float64: np.bincount returns int64 for an empty input even with weights, which
        # breaks the in-place accumulation across chunks (e.g. a subsample without events in a chunk)
        return np.bincount(idx, weights=weights, minlength=nb).astype(np.float64)

    out = {'n_raw': sums(), 'before': sums(w0)}
    for var, sf in sfs.items():
        out[f'after__{var}'] = sums(w0 * sf)
    return out


def compute_ratios(sums):
    """r per variation and jet bin for ONE (era, sample) from its (normalized) sums."""
    n_raw, before = sums['n_raw'], sums['before']
    out = {}
    for var in VARIATIONS:
        after = sums.get(f'after__{var}')
        if after is None:
            continue
        r_incl = before.sum() / after.sum() if after.sum() > 0 else 1.0
        ok = (n_raw >= MIN_EVENTS) & (after > 0) & (before > 0)
        r = np.where(ok, before / np.where(after > 0, after, 1.0), r_incl)
        out[var] = {'r': [float(v) for v in r], 'r_inclusive': float(r_incl),
                    'fallback_bins': [int(i) for i in np.flatnonzero(~ok)]}
    out['_n_raw'] = [int(n) for n in n_raw]
    return out


def _add(total, key, sums, norm):
    entry = total.setdefault(key, {})
    for name, arr in sums.items():
        arr = np.asarray(arr, dtype=np.float64) * (1.0 if name == 'n_raw' else norm)
        entry[name] = entry[name] + arr if name in entry else arr


def _sums_from_columns(cols):
    """Old outputs: rebuild the sums from the per-event columns of the calibration category."""
    data = {c: cols[f'events_{c}'].value for c in ['n_ak4jets', 'genWeight'] + OTHER_WEIGHTS
            + list(VARIATIONS.values()) if f'events_{c}' in cols}
    if 'n_ak4jets' not in data or 'btag_sf' not in data:
        return None
    w0 = np.sign(data['genWeight'])
    for col in OTHER_WEIGHTS:
        if col in data:
            w0 = w0 * data[col]
    sfs = {var: data[col] for var, col in VARIATIONS.items() if col in data}
    return calib_sums(data['n_ak4jets'], w0, sfs)


def read_calib(coffea_files, category='btag_calib', variations=None, samples=None):
    """{shape variation: {(era, sample): normalized sums}} from one or more .coffea outputs.
    variations=None reads all of them. Every (dataset, sample, variation) is taken once, so the
    nominal that is also present in the outputs of the jet-variation runs is not double counted."""
    from coffea.util import load
    from analysis_years import find_era
    from outputToPkl import get_norm

    total, seen = {}, set()

    def add(variation, dataset, sample, sums, genweight_df):
        if (variations and variation not in variations) or (dataset, sample, variation) in seen:
            return
        seen.add((dataset, sample, variation))
        era = find_era(dataset)
        _add(total.setdefault(variation, {}), (era, sample), sums,
             get_norm(sample, dataset, genweight_df, era))

    for coffea_file in coffea_files:
        filein = load(coffea_file)
        genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')

        if 'btag_calib' in filein:                       # accumulated by the processor
            for dataset, by_sample in filein['btag_calib'].items():
                for sample, by_var in by_sample.items():
                    if samples and sample not in samples:
                        continue
                    for variation, sums in by_var.items():
                        add(variation, dataset, sample, sums, genweight_df)
            continue

        print(f"  {coffea_file}: no 'btag_calib' accumulator, reading the '{category}' column category")
        for sample, datasets in filein['columns'].items():
            if samples and sample not in samples:
                continue
            for dataset, cats in datasets.items():
                for variation, cols in cats.get(category, {}).items():
                    sums = _sums_from_columns(cols)
                    if sums is None:
                        print(f"  [skip] {sample}/{dataset}: n_ak4jets or btag_sf missing in '{category}'")
                        continue
                    add(variation, dataset, sample, sums, genweight_df)
    return total


def build_table(sums_by_key, label='nominal'):
    table = {'njet_edges': NJET_EDGES}
    print(f"--- {label}")
    for (era, sample), sums in sorted(sums_by_key.items()):
        res = compute_ratios(sums)
        table.setdefault(era, {})[sample] = res
        nom = res.get('nominal', {})
        print(f"  {era:13} {sample:32} r(nominal) per jet bin: "
              + " ".join(f"{v:.3f}" for v in nom.get('r', []))
              + (f"   fallback bins {nom['fallback_bins']}" if nom.get('fallback_bins') else ""))
    return table


def build_tables(sums_by_variation):
    """Nominal table at the top level, the other shape variations under 'shape_variations'."""
    if 'nominal' not in sums_by_variation:
        raise ValueError("no nominal b-tag calibration sums in the inputs")
    table = build_table(sums_by_variation['nominal'])
    others = {v: build_table(s, v) for v, s in sorted(sums_by_variation.items()) if v != 'nominal'}
    if others:
        table['shape_variations'] = others
    return table


def apply_btag_norm(df, era, sample, table, variation='nominal'):
    """Multiply the b-tag SF columns of an analysis DataFrame by r(era, sample, njets).
    SF variations get their own r, so btag_sf<src>_up / btag_sf keeps the BTV-normalized ratio.
    For a shape variation the table of that variation is used (falls back to nominal)."""
    if variation != 'nominal':
        var_table = table.get('shape_variations', {}).get(variation)
        if var_table is None:
            print(f"  [WARNING] no b-tag normalization for shape variation {variation}; using nominal")
        else:
            table = var_table
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
    ap.add_argument('--input', '-i', required=True, nargs='+', help='.coffea output(s)')
    ap.add_argument('--output', '-o', default='btag_norm.json')
    ap.add_argument('--category', default='btag_calib',
                    help="column category to use for old outputs without the 'btag_calib' accumulator")
    ap.add_argument('--variations', nargs='+', default=None,
                    help='shape variations to compute r for (default: all found; nominal is always needed)')
    args = ap.parse_args()
    variations = None if args.variations is None else sorted(set(args.variations) | {'nominal'})
    table = build_tables(read_calib(args.input, args.category, variations))
    with open(args.output, 'w') as f:
        json.dump(table, f, indent=1)
    print(f"Wrote {args.output}")
