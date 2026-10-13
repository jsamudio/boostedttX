#                    #
##                  ##
######################
### Datacard shapes ##
######################
######################

import numpy as np
import pandas as pd
import json
import functools

from analysis_years import YEARS, PICKLE_DIR
from outvars import nominal_weight_factors

# Per-event correction factors applied on top of norm_weight * sign(genWeight).
# Columns missing for a process (e.g. topptWeight for non-ttbar samples) count as 1.
WEIGHT_COLUMNS = nominal_weight_factors


def getZhbbWeight(df_, year=None):
    """Total event weight. The year argument is kept for interface compatibility.
    FIX: the previous version used np.sign(genWeight * topptWeight), which applied only the
    sign of the top-pT weight (always +1), i.e. the top-pT reweighting was never applied."""
    w = df_['norm_weight'] * np.sign(df_['genWeight'])
    for col in WEIGHT_COLUMNS:
        if col in df_:
            w = w * df_[col]
    return w


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    """
    values = np.array(values, dtype='float')
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight, dtype='float')
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.asarray(np.cumsum(sample_weight) - 0.5 * sample_weight)
    if old_style:
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


class DataCardShapes():
    '''
    Handle format of datacard shapes.

    NN bin edges are signal-quantile based per reco pT channel. With common_binning=True
    (default) the edges are derived from the signal of all years together, so every year
    uses the same binning (simplest for combining years and for the EFT scaling).
    With common_binning=False each year gets its own edges.
    '''
    ref_samples = ['ttZ', 'ttH']
    quantiles = [0.0, .05, .25, .35, .50, .70, 1.0]

    def __init__(self, recopt_bins, recosdM_bins, n_NN_bins=10, nn='nnscore', isblind=True,
                 years=None, common_binning=True, file_dir=PICKLE_DIR, edges_file=None):
        # edges_file: JSON {year: {'Zhpt1': [...], ...}} from optimize_binning.py; replaces the
        # quantile edges (the first year in the file is used when a year is missing)
        self.fixed_edges = None
        if edges_file:
            with open(edges_file) as f:
                self.fixed_edges = json.load(f)
        self.years = list(years) if years is not None else list(YEARS)
        self.pt_bins = recopt_bins + [500]
        self.sdM_bins = recosdM_bins
        self.nn = nn
        self.isblind = isblind
        self.common_binning = common_binning
        self.file_dir = file_dir
        self.hist_dict = {}
        self.nn_edges = {}
        self.init_hist_funcs()

    def _signal_df(self, y):
        dfs = [pd.read_pickle(f"{self.file_dir.format(year=y)}/{s}.pkl") for s in self.ref_samples]
        df = pd.concat(dfs, axis='rows', ignore_index=True)
        df = df[(df[self.nn] >= 0.0) & (df['matchedGen_ZHbb_bb'] == True)].copy()
        df['ZH_pt'] = df['ZH_pt'].clip(self.pt_bins[0] + 1, self.pt_bins[-1] - 1)
        df['pt_bin'] = pd.cut(df['ZH_pt'], bins=self.pt_bins,
                              labels=[f'Zhpt{i}' for i in range(len(self.pt_bins) - 1)])
        df['_w'] = getZhbbWeight(df, y)
        return df

    def _edges(self, df, i):
        sub = df[df['pt_bin'] == f'Zhpt{i}']
        edges = weighted_quantile(sub[self.nn], self.quantiles, sub['_w'])
        edges[0], edges[-1] = 0, 1          # explicitly set bin edges to 0 and 1
        return edges

    def init_hist_funcs(self):
        if self.fixed_edges is not None:
            first = next(iter(self.fixed_edges.values()))
            for y in self.years:
                self.hist_dict[y], self.nn_edges[y] = {}, {}
                for i in range(1, len(self.pt_bins) - 1):
                    edges = np.asarray(self.fixed_edges.get(y, first)[f'Zhpt{i}'], dtype=float)
                    self.nn_edges[y][i] = edges
                    print(f'{y} Zhpt{i} NN edges (from file): {np.round(edges, 8).tolist()}')
                    mass = self.sdM_bins[i] if isinstance(self.sdM_bins, dict) else self.sdM_bins
                    self.hist_dict[y][i] = functools.partial(np.histogram2d, bins=[edges, mass])
            return
        frames = {y: self._signal_df(y) for y in self.years}
        pooled = pd.concat(list(frames.values()), ignore_index=True) if self.common_binning else None
        common = {}
        # pT bin 0 (0-200) is not used, start at index 1
        for y in self.years:
            df = pooled if self.common_binning else frames[y]
            self.hist_dict[y], self.nn_edges[y] = {}, {}
            for i in range(1, len(self.pt_bins) - 1):
                if self.common_binning and i in common:
                    edges = common[i]
                else:
                    edges = self._edges(df, i)
                    common[i] = edges
                self.nn_edges[y][i] = edges
                print(f'{y} Zhpt{i} NN edges: {np.round(edges, 8).tolist()}')
                mass = self.sdM_bins[i] if isinstance(self.sdM_bins, dict) else self.sdM_bins
                self.hist_dict[y][i] = functools.partial(np.histogram2d, bins=[edges, mass])

    def save_edges(self, path):
        """Write the NN bin edges per year and channel (e.g. for the EFT scaling's NN_BINS)."""
        out = {y: {f'Zhpt{i}': [float(v) for v in e] for i, e in d.items()}
               for y, d in self.nn_edges.items()}
        with open(path, 'w') as f:
            json.dump(out, f, indent=1)
        print(f'Saved NN bin edges to {path}')

    def __getitem__(self, y):  # building this like a dictionary
        try:
            return self.hist_dict[y]
        except KeyError:
            raise KeyError(f"{y} is not a valid top-level key!!!")