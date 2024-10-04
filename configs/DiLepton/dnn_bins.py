import pandas as pd
import numpy as np
import os
import sys
import awkward as ak
from glob import glob

from coffea.util import load
import argparse
import outvars

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Build datasets for NN training')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] >= 2) &
        (df_['ZH_bbvLscore'] >= 0.6) &
        (df_['n_ak4jets']   >= 5)             &
        #( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        #(df_['passNotHadLep'] == 1) & # might add
        (df_['ZH_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['ZH_M']        >= 50)            &
        (df_['ZH_M']        <= 200)
    )
    return base_cuts

pt_bins = [200, 300, 450, 500]

NN_vars = outvars.NN_vars
sig_vars = outvars.NN_vars+outvars.sig_vars
bkg_vars = outvars.NN_vars+outvars.bkg_vars

sig = ['ttHTobb__genMatch', 'ttHToNonbb__genMatch','TTZToBB__genMatch', 'TTZToQQ__genMatch', 'TTZToLLNuNu__genMatch']
bkg = ["TTbb_Hadronic__tt+B",
       "TTbb_SemiLeptonic__tt+B",
       "TTbb_2L2Nu__tt+B",
       "TTToSemiLeptonic__tt+LF",
       "TTToSemiLeptonic__tt+C",
       "TTTo2L2Nu__tt+LF",
       "TTTo2L2Nu__tt+C",
       "TTToHadronic__tt+LF",
       "TTToHadronic__tt+C"]

genmatchreq = 'matchedGen_ZHbb_bb'

def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.asarray(np.cumsum(sample_weight) - 0.5 * sample_weight)
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)
#

class DNN_bins:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    #bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore', 'newgenm_NN','norm_weight', 'genWeight']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df])
        self.nn_bins = self.get_NN_bins()


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+[genmatchreq]:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+[genmatchreq])
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make signal df

        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+['tt_type']:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+['tt_type'])
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        s_df = s_df[s_df[genmatchreq] == True]
        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

        return s_df, b_df

    def get_NN_bins(self):
        cuts = dnn_cut(self.sb_df)
        genm = (self.sb_df[genmatchreq] == True)
        self.sb_df['ZH_pt'].clip(pt_bins[0]+1,pt_bins[-1]-1, inplace=True)
        # currently not binning by pt...
        quantiles = [0.0, .05, .25, .35, .50, .70, 1.0]
        for i in range(len(pt_bins)-1):
            nn_df = np.asarray(self.sb_df[ cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) &
                (self.sb_df['ZH_pt'] > pt_bins[i]) & (self.sb_df['ZH_pt'] < pt_bins[i+1])]['newgenm_NN'].to_numpy(), dtype = float)
            norm_weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) &
                (self.sb_df['ZH_pt'] > pt_bins[i]) & (self.sb_df['ZH_pt'] < pt_bins[i+1])]['norm_weight'].to_numpy(), dtype = float)
            weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) &
                (self.sb_df['ZH_pt'] > pt_bins[i]) & (self.sb_df['ZH_pt'] < pt_bins[i+1])]['genWeight'].to_numpy(), dtype = float)
            norm_weight = norm_weight * np.sign(weight) / 1000

            nn_bins = weighted_quantile(nn_df, quantiles, norm_weight)
            nn_bins[0], nn_bins[-1] = 0,1
            print(nn_bins)
        return nn_bins




if __name__ == '__main__':
    _ = DNN_bins()
