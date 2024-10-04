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
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
from plotUtils import get_sigbkg

parser = argparse.ArgumentParser(description='Build datasets for NN training')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] >= 2) &
        (df_['ZH_bbvLscore'] >= 0.6) &
        (df_['n_ak4jets']   >= 5)             &
        #FIXME these need to be added in

        #( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        #(df_['passNotHadLep'] == 1) & # might add
        (df_['ZH_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['ZH_M']        >= 50)            &
        (df_['ZH_M']        <= 200)
    )
    return base_cuts

def plot_cut(df_):
    base_cuts = (
            (df_['ZH_pt'] > 300)     &
            (df_['ZH_M'] > 75)       &
            (df_['ZH_M'] < 145)      &
            (df_['newgenm_NN'] <= 1) &
            (df_['newgenm_NN'] > 0.)
    )
    return base_cuts



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

class DNN_datasets:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    #bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore', 'newgenm_NN','norm_weight', 'genWeight', 'topptWeight']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = get_sigbkg(filein, sig, bkg, self.dnn_vars+self.cut_vars, genmatchreq)
        print(self.s_df)
        self.s_df = self.s_df[(genmatchreq)]
        print(self.b_df)
        #self.sb_df = pd.concat([self.s_df.reset_index(drop=True),self.b_df.reset_index(drop=True)], ignore_index=True)
        self.sb_df = pd.concat([self.s_df,self.b_df], ignore_index=True)
        print(self.sb_df)
        #self.sb_df = self.s_df.append(self.b_df)
        self.nn_bins = self.get_NN_bins()
        self.plot_dnnHist()

    def plot_dnnHist(self):
        fig, (ax, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        print(np.unique(self.sb_df.process, return_counts=True))
        cuts = dnn_cut(self.sb_df)
        plotcut = plot_cut(self.sb_df)
        sigs = ['ttZ', 'ttH', 'old_ttZbb']
        sumS = []
        sumB = []

        for i in self.sb_df.process.unique():
            norm_weight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
            weight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            topptWeight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)

            norm_weight = (topptWeight * norm_weight * np.sign(weight))

            n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
                    histtype='step', range= (0,1), label=f'{i}', weights=norm_weight)

            #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
            #        histtype='step', range= (0,1), label=f'{i}')
            if i in sigs:
                sumS.append(n)
            else:
                sumB.append(n)
            print(i, np.sum(n))

        #vals = pd.DataFrame(np.concatenate(sumS))
        #print(vals)
        print(np.sum(sumS,axis=0), np.sum(sumB,axis=0))
        sumS = np.sum(sumS,axis=0)
        sumB = np.sum(sumB,axis=0)
        bin_c = (bins[1:]+bins[:-1])/2
        ax2.errorbar(x=bin_c, y = sumS/np.sqrt(sumB), xerr=(bins[1:]-bins[:-1])/2,
                fmt='.', color='k', label=r'S/$\sqrt{\mathrm{B}}$')
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='both', direction='in', top=True, right=True)
        ax2.yaxis.set_label_coords(-0.07,0.35)
        ax2.set_ylabel(r'$\mathrm{S/}\sqrt{\mathrm{B}}$')
        ax2.grid(True)


        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.set_yscale('log')
        ax.set_xlim([bins[0],bins[-1]])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=8)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0,
            wspace=0.2)

        plt.savefig("testDNN.pdf")

    def get_NN_bins(self):
        cuts = dnn_cut(self.sb_df)
        genm = (self.sb_df[genmatchreq] == True)
        # currently not binning by pt...
        nn_df = np.asarray(self.sb_df[ cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) & (self.sb_df['ZH_pt'] > 300)]['newgenm_NN'].to_numpy(), dtype = float)
        norm_weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) & (self.sb_df['ZH_pt'] > 300)]['norm_weight'].to_numpy(), dtype = float)
        weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) & (self.sb_df['ZH_pt'] > 300)]['genWeight'].to_numpy(), dtype = float)
        topptWeight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['newgenm_NN'] <= 1) & (self.sb_df['newgenm_NN'] > 0.) & (self.sb_df['ZH_pt'] > 300)]['topptWeight'].to_numpy(), dtype = float)
        norm_weight = (topptWeight * norm_weight * np.sign(weight))
        quantiles = [0.0, .05, .25, .35, .50, .70, 1.0]

        nn_bins = weighted_quantile(nn_df, quantiles, norm_weight)
        nn_bins[0], nn_bins[-1] = 0,1
        print(nn_bins)
        nn_bins = [0., 0.083, 0.431, 0.560, 0.736, 0.865, 1. ]
        return nn_bins




if __name__ == '__main__':
    _ = DNN_datasets()
