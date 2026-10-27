"""
Convert a PocketCoffea .coffea output into per-process pickles, one directory per analysis year:
    pickled/Nominal<year>/<process>.pkl      (see PICKLE_DIR in analysis_years.py)

The era of each dataset is read from its name (e.g. TTToHadronic_2016_PreVFP) and the
normalization uses that era's luminosity. Eras are mapped to analysis years via ERA_TO_YEAR
(by default 2016_PreVFP + 2016_PostVFP -> '2016'); every event keeps its era in the 'era' column.

Usage:
    python outputToPkl.py -i output_signal.coffea output_bkg.coffea output_minor.coffea \
                          [--btag-norm btag_norm.json] [--variations nominal|all|<name> ...]

With --btag-norm, the b-tag SF columns are multiplied by the BTV normalization factor
r(era, sample, n_ak4jets) from btag_norm.py.

Shape variations (e.g. a PC_JEC_VARIATION run) are exported as <process>_<variation>.pkl,
e.g. ttZ_AK4AK8_JES_AbsoluteStatUp.pkl, which is the naming makeDatacard.jet_systematics expects.
"""
import gc
import json
import argparse
import os
import re

import numpy as np
import pandas as pd
from coffea.util import load

import outvars
from analysis_years import ERA_RE, find_era, ERA_TO_YEAR, LUMI_FB, PICKLE_DIR
from samples import PICKLE_SIGNAL, PICKLE_BACKGROUND
from spanet_models import ACTIVE_MODELS, score_columns

ERA_SUFFIX = re.compile(rf'_{ERA_RE}$')

SIG = PICKLE_SIGNAL          # (sub)samples exported, see samples.py
BKG = PICKLE_BACKGROUND
CATEGORY = 'btag_mask'

PROCESSES = ['ttZ', 'ttH', 'TTBar', 'tt_B', 'QCD', 'VJets', 'SingleTop', 'ttX']

# Process-label renaming after reading (background samples -> datacard process names)
RELABEL = {'WJetsToLNu': 'VJets', 'DYJetsToLL': 'VJets'}

# 13 TeV cross sections in pb (Run 2). Run 3 years need 13.6 TeV values.
XSEC_DICT = {
    "TTZToLLNuNu": 0.2529, "ttHTobb": 0.2934, "ttHToNonbb": 0.215,
    "TTZToQQ": 0.5297, "TTZToBB": 0.1157, "TTToHadronic": 380.095,
    "TTToSemiLeptonic": 364.018, "TTTo2L2Nu": 88.29, "TTbb_Hadronic": 9.286,
    "TTbb_SemiLeptonic": 8.168, "TTbb_2L2Nu": 1.807, "TTWJetsToLNu": 0.2163,
    "TTWJetsToQQ": 0.4432, "TTTT": 0.008091, "TTGJets": 3.774,
    "ST_tW_top": 39.65, "ST_tW_antitop": 39.65, "ST_s_channel": 3.608,
    "ST_t_channel_top": 134.2, "ST_t_channel_antitop": 80.0,
    "QCD_HT500to700": 29980.0, "QCD_HT700to1000": 6334.0, "QCD_HT1000to1500": 1088.0,
    "QCD_HT1500to2000": 99.11, "QCD_HT2000toInf": 20.23,
    "WJetsToLNu_HT-70To100": 1524.0, "WJetsToLNu_HT-100To200": 1256.0,
    "WJetsToLNu_HT-200To400": 335.5, "WJetsToLNu_HT-400To600": 45.25,
    "WJetsToLNu_HT-600To800": 10.97, "WJetsToLNu_HT-800To1200": 4.933,
    "WJetsToLNu_HT-1200To2500": 1.074, "WJetsToLNu_HT-2500ToInf": 0.03216,
    "DYJetsToLL_HT-70to100": 175.3, "DYJetsToLL_HT-100to200": 139.4,
    "DYJetsToLL_HT-200to400": 38.4, "DYJetsToLL_HT-400to600": 5.174,
    "DYJetsToLL_HT-600to800": 1.258, "DYJetsToLL_HT-800to1200": 0.5765,
    "DYJetsToLL_HT-1200to2500": 0.1331, "DYJetsToLL_HT-2500toInf": 0.002987
}

MODELS = ACTIVE_MODELS       # see spanet_models.py


def get_norm(sample, dataset, genweight_df, era):
    """Per-event normalization: lumi(era) * xsec / sum of sign(genWeight). Fails loudly."""
    name = ERA_SUFFIX.sub('', dataset)
    key = name if name in XSEC_DICT else sample.split('__')[0]
    if key not in XSEC_DICT:
        raise KeyError(f"No cross section for dataset '{dataset}' (sample '{sample}')")
    if dataset not in genweight_df.index:
        raise KeyError(f"No sum_signOf_genweights for dataset '{dataset}'")
    if era not in LUMI_FB:
        raise KeyError(f"No luminosity configured for {era} (LUMI_FB in analysis_years.py)")
    return LUMI_FB[era] * 1000.0 * XSEC_DICT[key] / genweight_df.loc[dataset, 0]


class OutputToPkl:
    """The input files are read one at a time: only the exported columns of each file are kept
    (as DataFrames, per variation) and the file is freed before the next one is opened, so the
    peak memory is the largest single .coffea file plus the exported columns."""

    def __init__(self, input_files, btag_norm=None):
        if isinstance(input_files, str):
            input_files = [input_files]
        self.input_files = list(input_files)
        self.btag_norm = json.load(open(btag_norm)) if btag_norm else None

        spanet_vars = [c for m in MODELS for c in score_columns(m)]
        # PDF members exist only for signal and tt+jets (outvars.pdf_vars): read when present
        self.sig_vars = (outvars.common_vars + outvars.sig_vars + outvars.weight_vars + spanet_vars + ["xbb_sf"]
                         + outvars.pdf_vars)
        self.bkg_vars = (outvars.common_vars + outvars.weight_vars + spanet_vars + ["xbb_sf"] + outvars.pdf_vars
                         + outvars.bkg_vars)

    @staticmethod
    def variations_in(filein):
        """Shape variations present in one input (nominal first)."""
        found = set()
        for datasets in filein['columns'].values():
            for cats in datasets.values():
                found.update(cats.get(CATEGORY, {}).keys())
        return sorted(found, key=lambda v: (v != 'nominal', v))

    def read_samples(self, filein, genweight_df, samples, variables, variation='nominal'):
        """Read every dataset of the given samples in one input, with per-dataset year and normalization."""
        frames = []
        for sample in samples:
            if sample not in filein['columns']:
                continue
            for dataset in filein['columns'][sample].keys():
                by_var = filein['columns'][sample][dataset].get(CATEGORY, {})
                if variation not in by_var:
                    continue
                cols = by_var[variation]
                era = find_era(dataset)
                year = ERA_TO_YEAR[era]
                norm = get_norm(sample, dataset, genweight_df, era)

                tmp, missing = {}, []
                for var in variables:
                    if var == 'norm_weight':
                        continue
                    key = f'events_{var}'
                    if key in cols:
                        tmp[var] = cols[key].value
                    elif var not in outvars.pdf_vars + outvars.bkg_vars:   # these exist only for some samples
                        missing.append(var)
                if missing:
                    print(f"  [WARNING] {sample}/{dataset}: missing {missing}")
                if tmp:
                    df = pd.DataFrame(tmp)
                    if self.btag_norm is not None:
                        from btag_norm import apply_btag_norm
                        if 'n_ak4jets' in df:
                            df = apply_btag_norm(df, era, sample, self.btag_norm, variation)
                        else:
                            print(f"  [WARNING] {sample}/{dataset}: no n_ak4jets, b-tag normalization skipped")
                    df['norm_weight'] = norm
                    df['era'] = era
                    df['year'] = year
                    frames.append(df)
                    print(f"  {sample}/{dataset} [{variation}]: era={era} -> year={year}, "
                          f"{len(df)} events, norm={norm:.4g}")
        return frames

    def run(self, variations=('nominal',), skip_nominal=False):
        want_all = 'all' in variations
        frames = {}              # variation -> {'sig': [DataFrame], 'bkg': [DataFrame]}
        for path in self.input_files:
            print(f"Reading {path}")
            filein = load(path)
            genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')
            todo = self.variations_in(filein) if want_all else list(variations)
            if skip_nominal:
                todo = [v for v in todo if v != 'nominal']
            for variation in todo:
                f = frames.setdefault(variation, {'sig': [], 'bkg': []})
                f['sig'] += self.read_samples(filein, genweight_df, SIG, self.sig_vars, variation)
                f['bkg'] += self.read_samples(filein, genweight_df, BKG, self.bkg_vars, variation)
            del filein, genweight_df
            gc.collect()

        for variation in sorted(frames, key=lambda v: (v != 'nominal', v)):
            print('-' * 60)
            print(f"Variation: {variation}")
            f = frames.pop(variation)
            s_df = pd.concat(f['sig'], ignore_index=True) if f['sig'] else pd.DataFrame()
            b_df = pd.concat(f['bkg'], ignore_index=True) if f['bkg'] else pd.DataFrame()
            del f
            self.run_variation(variation, s_df, b_df)

    def run_variation(self, variation, s_df, b_df):
        if not b_df.empty:
            b_df['process'] = b_df['process'].replace(RELABEL)
            b_df = b_df[b_df['process'].isin(['TTBar', 'tt_B', 'QCD', 'VJets', 'ttX', 'SingleTop'])]

        # Placeholder data (copy of the signal) until real data is added; nominal only
        data_df = s_df.copy() if (variation == 'nominal' and not s_df.empty) else pd.DataFrame()
        if not data_df.empty:
            data_df['process'] = 'data_obs'
            data_df['norm_weight'] = 1.0

        sb_df = pd.concat([s_df, b_df, data_df], ignore_index=True)
        if sb_df.empty:
            print(f"  nothing to export for {variation}")
            return

        for m in MODELS:
            sig, ttbb, ttlf = sb_df[f"sig_score_{m}"], sb_df[f"ttbb_score_{m}"], sb_df[f"ttlf_score_{m}"]
            denom = sig + ttbb + ttlf
            sb_df[f"disc_{m}"] = np.where(denom > 0, sig / np.where(denom > 0, denom, 1.0), 0.0)

        self.export(sb_df, suffix='' if variation == 'nominal' else f'_{variation}')

    @staticmethod
    def export(sb_df, suffix=''):
        for year, ydf in sb_df.groupby('year'):
            out_dir = PICKLE_DIR.format(year=year)
            os.makedirs(out_dir, exist_ok=True)
            for proc in PROCESSES + ['data_obs']:
                export_df = ydf[ydf['process'] == proc]
                if not export_df.empty:
                    path = os.path.join(out_dir, f"{proc}{suffix}.pkl")
                    export_df.to_pickle(path)
                    print(f"Exported: {path} ({len(export_df)} events)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Change .coffea output(s) to per-year .pkl files')
    parser.add_argument('--input', '-i', type=str, nargs='+', required=True, help='Input .coffea file(s)')
    parser.add_argument('--btag-norm', type=str, default=None,
                        help='JSON from btag_norm.py; applies the b-tag shape SF normalization')
    parser.add_argument('--variations', nargs='+', default=['nominal'],
                        help="shape variations to export: 'nominal' (default), 'all', or explicit names")
    parser.add_argument('--skip-nominal', action='store_true',
                        help='do not (re)export the nominal pickles, e.g. for the output of a JEC-variation run')
    args = parser.parse_args()
    OutputToPkl(args.input, btag_norm=args.btag_norm).run(args.variations, skip_nominal=args.skip_nominal)
