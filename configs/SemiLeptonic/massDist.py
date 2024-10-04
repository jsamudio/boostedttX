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

parser = argparse.ArgumentParser(description='Generate mass histogram')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def dnn_cut(df_):
    base_cuts = (
        (df_['n_ak4jets']   >= 5)       &
        #(df_['ZH_pt']       >= 200)     &
        (df_['ZH_pt']       >  450)     &
        #(df_['newgenm_NN']       >  0.7)
        (df_['newgenm_NN']       >  0.94)
    )
    return base_cuts


NN_vars = outvars.NN_vars
sig_vars = outvars.NN_vars+outvars.sig_vars
bkg_vars = outvars.NN_vars+outvars.bkg_vars

sig = ['ttHTobb__genMatch',
       'ttHTobb__non_genMatch',
       'ttHToNonbb__genMatch',
       'ttHToNonbb__non_genMatch',
       'TTZToBB__genMatch',
       'TTZToBB__non_genMatch',
       'TTZToQQ__genMatch',
       'TTZToQQ__non_genMatch',
       'TTZToLLNuNu__genMatch',
       'TTZToLLNuNu__non_genMatch']

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

class DNN_datasets:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    #bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore', 'newgenm_NN','norm_weight', 'genWeight', 'topptWeight']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df])
        self.plot_dnnHist()


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]

        genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')
        print(genweight_df[0]['TTToSemiLeptonic__2017'])

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+[genmatchreq]:
                    if (var == 'norm_weight'):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+[genmatchreq])
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make bkg df

        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+['tt_type']:
                    if ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+['tt_type'])
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        #s_df = s_df[s_df[genmatchreq] == True]
        s_df = s_df.replace('old_ttZbb', 'ttZ')
        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

        return s_df, b_df

    def plot_dnnHist(self):
        fig, ax = plt.subplots()
        print(np.unique(self.sb_df.process, return_counts=True))
        cuts = dnn_cut(self.sb_df)
        sigs = ['ttZ', 'ttH']
        bkgs = ['tt_B', 'TTBar']
        sumS = []
        sumB = []

        for i in bkgs+sigs:
            norm_weight = ((np.asarray(self.sb_df[cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float))/(41.529))*137.596
            norm_weight2 = ((np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float))/(41.529))*137.596
            weight = np.asarray(self.sb_df[cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            weight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            topptWeight = np.asarray(self.sb_df[cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
            topptWeight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)

            norm_weight = (topptWeight * norm_weight * np.sign(weight))
            norm_weight2 = (topptWeight2 * norm_weight2 * np.sign(weight2))

            if i in bkgs:
                n, bins, patches = ax.hist(self.sb_df['ZH_M'][cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5), stacked=True,
                    histtype='stepfilled', label=f'{i}', weights=norm_weight)
                print("okay")
            else:
                n, bins, patches = ax.hist(self.sb_df['ZH_M'][cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5), stacked=False,
                    histtype='step', label=f'{i} x 5', weights=norm_weight*5)
                n, bins, patches = ax.hist(self.sb_df['ZH_M'][(self.sb_df[genmatchreq] == True) & cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5), stacked=False,
                    histtype='step', label=f'{i} GenMatched x 5', weights=norm_weight2*5)
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
        #ax2.errorbar(x=bin_c, y = sumS/np.sqrt(sumB), xerr=(bins[1:]-bins[:-1])/2,
        #        fmt='.', color='k', label=r'S/$\sqrt{\mathrm{B}}$')
        #ax2.xaxis.set_minor_locator(AutoMinorLocator())
        #ax2.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        #ax2.yaxis.set_minor_locator(AutoMinorLocator())
        #ax2.tick_params(which='both', direction='in', top=True, right=True)
        #ax2.yaxis.set_label_coords(-0.07,0.35)
        #ax2.set_ylabel(r'$\mathrm{S/}\sqrt{\mathrm{B}}$')
        #ax2.grid(True)
        tex_x_corr = 0.55
        #fig.text(tex_x_corr,0.59, r'200 < p_T < 300 GeV', fontsize=7)
        #fig.text(tex_x_corr,0.59, r'300 < p_T < 450 GeV', fontsize=7)
        fig.text(tex_x_corr,0.59, r'p_T > 450 GeV', fontsize=7)
        fig.text(tex_x_corr,0.53, r'DNN score > 0.94', fontsize=7)

        ax.set_ylabel("Events / 5 GeV")
        ax.set_xlabel(r"$m_{PNet}^\text{Z/H cand.}$")

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='both', direction='in', top=True, right=True)
        #ax.set_yscale('log')
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
        ax.set_ylim(0, 45)

        plt.savefig("pdf/mass450.pdf")


if __name__ == '__main__':
    _ = DNN_datasets()
