"""
Run 3 copy of Run2/outputToPkl.py: convert PocketCoffea .coffea outputs into per-process pickles,
one directory per analysis year:
    pickled/Nominal<year>/<process>.pkl      (see PICKLE_DIR in analysis_years.py)

The era of each dataset is read from its name (e.g. TTtoLNu2Q_2024) and the normalization uses that
era's luminosity; every event keeps its era in the 'era' column.

Usage:
    python outputToPkl.py -i outputs/nominal/signal/output_all.coffea outputs/nominal/bkg/output_all.coffea \
                          outputs/nominal/minor/output_all.coffea [--variations nominal|all|<name> ...]

Differences from Run 2:
  * 13.6 TeV cross sections (XSEC_DICT, by dataset name), Run 3 lumi (analysis_years.LUMI_FB)
  * no --btag-norm: the fixed-WP b-tag SFs (method 1a, with our own efficiency maps) need no
    normalization factor, unlike the Run 2 shape SFs (btag_norm.py)
  * no tH SM reweighting (th_reweighting.py): the 2024 THQ / THW samples are taken as SM
  * VV is its own process (no diboson sample in Run 2)

Shape variations (e.g. a PC_JEC_VARIATION run) are exported as <process>_<variation>.pkl,
e.g. ttZ_AK4PFPuppi_JES_TotalUp.pkl, which is the naming makeDatacard.jet_systematics expects.
"""
import gc
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

TTBB_SCALE_FIX_PREFIX = "TTbb_"      # = custom_scale_factors.TTBB_SCALE_FIX_PREFIX (no coffea import here)


def is_ttbb_scale_fix(sample):
    return sample.startswith(TTBB_SCALE_FIX_PREFIX)

ERA_SUFFIX = re.compile(rf'_{ERA_RE}$')
# ttbar / ttbb modelling-variation samples (samples.TTBAR_VARIATION_SAMPLES): nominal cross section
VARIATION_SUFFIX = re.compile(r'_(hdampUP|hdampDOWN|TuneCP5up|TuneCP5down)$')

SIG = PICKLE_SIGNAL          # (sub)samples exported, see samples.py
BKG = PICKLE_BACKGROUND
CATEGORY = 'btag_mask'

PROCESSES = ['ttZ', 'ttH', 'TTBar', 'tt_B', 'QCD', 'VJets', 'VV', 'SingleTop', 'ttX']

# Process-label renaming after reading (samples.PROCESS_LABELS already gives the datacard names)
RELABEL = {}

# 13.6 TeV cross sections in pb, by dataset name (without the era suffix). PROVISIONAL: the values of
# the dataset definitions (datasets/*.json metadata, the ones the processor's norm_weight uses);
# their sources are still to be documented / replaced by the reference values (see README, Run 3
# datacards). Run 2 for comparison (13 TeV): TTToSemiLeptonic 364.018, TTbb_SemiLeptonic 8.168.
XSEC_DICT = {
    # signal
    "TTH-Hto2B": 0.331, "TTH-HtoNon2B": 0.238,
    "TTZ-ZtoQQ-1Jets": 0.7865, "TTLL_Bin-MLL-4to50": 0.05425, "TTLL_Bin-MLL-50": 0.1188, "TTNuNu": 0.2250,
    # ttbar (5FS) and ttbb (4FS)
    "TTto4Q": 419.31, "TTtoLNu2Q": 405.4604, "TTto2L2Nu": 97.9016,
    "TTBBto4Q": 19.15297, "TTBBtoLNu2Q": 17.36875, "TTBBto2L2Nu": 3.74766,
    # ttX
    "TTWJets": 0.4678, "TTTT": 0.009652, "THW": 0.0172, "THQ": 0.08362, "TZQ_ll": 0.07968,
    # single top
    "ST_tW_top-4Q": 19.9533, "ST_tW_top-LNu2Q": 19.29405, "ST_tW_top-2L2Nu": 4.6587,
    "ST_tW_antitop-4Q": 19.9533, "ST_tW_antitop-LNu2Q": 19.29405, "ST_tW_antitop-2L2Nu": 4.6587,
    "ST_top_s_lep": 2.360, "ST_antitop_s_lep": 1.477,
    "ST_top_t_2Q": 97.74, "ST_top_t_LNu": 47.26, "ST_antitop_t_2Q": 58.78, "ST_antitop_t_LNu": 28.42,
    # V+jets, VV
    "WtoLNu-4Jets_Bin-1J": 9141.0, "WtoLNu-4Jets_Bin-2J": 2931.0, "WtoLNu-4Jets_Bin-3J": 864.6,
    "WtoLNu-4Jets_Bin-4J": 417.8,
    "DYto2E_MLL-10to50": 5743.0, "DYto2E_MLL-50": 1827.0, "DYto2Mu_MLL-10to50": 5803.0,
    "DYto2Mu_MLL-50": 1815.0, "DYto2Tau_MLL-10to50": 5806.0, "DYto2Tau_MLL-50": 1816.0,
    "WW": 125.8, "WZ": 54.7, "ZZ": 16.7,
    # QCD
    "QCD-HT_40_70": 312300000.0, "QCD-HT_70_100": 58470000.0, "QCD-HT_100_200": 25310000.0,
    "QCD-HT_200_400": 1960000.0, "QCD-HT_400_600": 97400.0, "QCD-HT_600_800": 13560.0,
    "QCD-HT_800_1000": 3010.0, "QCD-HT_1000_1200": 890.3, "QCD-HT_1200_1500": 384.8,
    "QCD-HT_1500_2000": 127.3, "QCD-HT_2000": 26.26,
}
# ttbar / ttbb modelling-variation samples (samples.TTBAR_VARIATION_SAMPLES, datasets named
# <sample>_<variation>_<era>): the cross section of the nominal dataset of the sample
VARIATION_NOMINAL = {"TTToHadronic": "TTto4Q", "TTToSemiLeptonic": "TTtoLNu2Q", "TTTo2L2Nu": "TTto2L2Nu",
                     "TTbb_Hadronic": "TTBBto4Q", "TTbb_SemiLeptonic": "TTBBtoLNu2Q", "TTbb_2L2Nu": "TTBBto2L2Nu"}


def xsec_key(sample, dataset):
    """XSEC_DICT key of a dataset: its name without the era, or for a variation sample the nominal dataset."""
    name = ERA_SUFFIX.sub('', dataset)
    if name in XSEC_DICT:
        return name
    base = VARIATION_SUFFIX.sub('', name)
    nominal = VARIATION_NOMINAL.get(base, VARIATION_NOMINAL.get(VARIATION_SUFFIX.sub('', sample.split('__')[0])))
    if nominal is None or nominal not in XSEC_DICT:
        raise KeyError(f"No cross section for dataset '{dataset}' (sample '{sample}')")
    return nominal


MODELS = ACTIVE_MODELS       # see spanet_models.py


def get_norm(sample, dataset, genweight_df, era, theory_sums=None):
    """Per-event normalization: lumi(era) * xsec / sum of sign(genWeight). Fails loudly.
    4FS ttbb (nominal moved to (mu0, mu0), custom_scale_factors.ttbb_scale_ref): the denominator is
    sum sign(genWeight) * ttbb_scale_ref over all generated events (theory_sums[sample][dataset]['all']
    ['nominal']), so the yield keeps the cross section of XSEC_DICT (as the tH SM reweighting of Run 2)."""
    key = xsec_key(sample, dataset)
    if era not in LUMI_FB:
        raise KeyError(f"No luminosity configured for {era} (LUMI_FB in analysis_years.py)")
    base = sample.split('__')[0]
    if is_ttbb_scale_fix(base):
        try:
            denominator = theory_sums[base][dataset]['all']['nominal']
        except (KeyError, TypeError):
            raise KeyError(f"No theory_sums for ttbb dataset '{dataset}' (sample '{base}'): needed for the "
                           f"(mu0, mu0) normalization (custom_scale_factors.ttbb_scale_ref)") from None
        return LUMI_FB[era] * 1000.0 * XSEC_DICT[key] / denominator
    if dataset not in genweight_df.index:
        raise KeyError(f"No sum_signOf_genweights for dataset '{dataset}'")
    return LUMI_FB[era] * 1000.0 * XSEC_DICT[key] / genweight_df.loc[dataset, 0]


class OutputToPkl:
    """The input files are read one at a time: only the exported columns of each file are kept
    (as DataFrames, per variation) and the file is freed before the next one is opened, so the
    peak memory is the largest single .coffea file plus the exported columns."""

    def __init__(self, input_files):
        if isinstance(input_files, str):
            input_files = [input_files]
        self.input_files = list(input_files)

        spanet_vars = [c for m in MODELS for c in score_columns(m)]
        # PDF members exist only for signal and tt+jets (outvars.pdf_vars): read when present
        self.sig_vars = (outvars.common_vars + outvars.sig_vars + outvars.weight_vars + spanet_vars
                         + outvars.pdf_vars + outvars.event_id_vars)
        self.bkg_vars = (outvars.common_vars + outvars.weight_vars + spanet_vars + outvars.pdf_vars
                         + outvars.bkg_vars + outvars.event_id_vars)

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
                norm = get_norm(sample, dataset, genweight_df, era, filein.get('theory_sums'))

                tmp, missing = {}, []
                for var in variables:
                    if var == 'norm_weight':
                        continue
                    key = f'events_{var}'
                    if key in cols:
                        tmp[var] = cols[key].value
                    elif var not in outvars.pdf_vars + outvars.bkg_vars + outvars.event_id_vars:
                        # pdf/bkg columns exist only for some samples; event ids not in older outputs
                        missing.append(var)
                if missing:
                    print(f"  [WARNING] {sample}/{dataset}: missing {missing}")
                if tmp:
                    df = pd.DataFrame(tmp)
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
            b_df = b_df[b_df['process'].isin(['TTBar', 'tt_B', 'QCD', 'VJets', 'VV', 'ttX', 'SingleTop'])]

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
    parser.add_argument('--variations', nargs='+', default=['nominal'],
                        help="shape variations to export: 'nominal' (default), 'all', or explicit names")
    parser.add_argument('--skip-nominal', action='store_true',
                        help='do not (re)export the nominal pickles, e.g. for the output of a JEC-variation run')
    args = parser.parse_args()
    OutputToPkl(args.input).run(args.variations, skip_nominal=args.skip_nominal)
