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

import mplhep as hep
hep.style.use("CMS")

parser = argparse.ArgumentParser(description='Generate mass histogram')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')
parser.add_argument('--pt', type=str, help='l, m, h')

args=parser.parse_args()

filein = load(args.input)
ptrange = args.pt

dirname = os.path.dirname(args.input)

if ptrange == 'l':
    ptlow = 200
    pthigh = 300
if ptrange == 'm':
    ptlow = 300
    pthigh = 450
if ptrange == 'h':
    ptlow = 450
    pthigh = np.inf

def dnn_cut(df_):
    base_cuts = (
        (df_['n_ak4jets']   >= 5)       &
        (df_['n_b_outZH'] == 2) &
        (df_['ZH_pt']       >= ptlow)     &
        (df_['ZH_pt']       <  pthigh) &
        (df_['ZH_bbvLscore'] >= 0.9105) &
        (df_['ZH_M'] >= 50) &
        (df_['ZH_M'] <= 200)
        #(df_['newgenm_NN']       >  0.7)
        #(df_['sigVsbkg']       >  0.8758)
        #(df_['sigVsbkg']       <  0.8758)
        #(df_['tthbb']       >  0.76)
    )
    return base_cuts


NN_vars = outvars.NN_vars
sig_vars = outvars.NN_vars+outvars.sig_vars
bkg_vars = outvars.NN_vars+outvars.bkg_vars

sig = ['ttHTobb__genMatch',
       'ttHTobb__non_genMatch',
       #'ttHToNonbb__genMatch',
       #'ttHToNonbb__non_genMatch',
       'TTZToBB__genMatch',
       'TTZToBB__non_genMatch',
       #'TTZToQQ__genMatch',
       #'TTZToQQ__non_genMatch',
       #'TTZToLLNuNu__genMatch',
       #'TTZToLLNuNu__non_genMatch'
      ]


bkg = ["TTbb_Hadronic__tt+B",
       "TTbb_SemiLeptonic__tt+B",
       "TTbb_2L2Nu__tt+B",
       #"TTToSemiLeptonic__tt+LF",
       #"TTToSemiLeptonic__tt+C",
       #"TTTo2L2Nu__tt+LF",
       #"TTTo2L2Nu__tt+C",
       #"TTToHadronic__tt+LF",
       #"TTToHadronic__tt+C"
       ]
vjets = [#"WJetsToLNu_HT",
       #"DYJetsToLL_HT",
       ]

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

class DNN_datasets:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    #bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore', 'norm_weight', 'genWeight', 'topptWeight','ttzbb', 'ttbb', 'ttlf', 'tthbb']# ,'signal', 'ZH_sdm'] 'min_wpart_ZH_dR', 'max_wpart_ZH_dR', 'n_wpart_ZH_dR_0p4', 'n_wpart_ZH_dR_0p8', 'n_wpart_ZH_dR_1p2', 'min_topb_ZH_dR', 'max_topb_ZH_dR', 'n_topb_ZH_dR_0p4', 'n_topb_ZH_dR_0p8', 'n_topb_ZH_dR_1p2']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df])
        self.nn_bins = self.get_NN_bins()
        self.plot_dnnHist(self.nn_bins)


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]

        genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')
        #print(genweight_df[0]['TTToSemiLeptonic__2017'])

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            #print(filein['columns'].keys())
            for j in filein['columns'][f'{i}'].keys():
                
                for var in pre_vars+[genmatchreq]:
                    if (var == 'norm_weight'):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb','tthbb', 'ttbb', 'ttlf', 'ttcc','signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputH_{var}'].value.tolist()
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
                #print(filein['columns'][f'{i}'][f'{j}']['btag_mask'].keys())
                for var in pre_vars+['ttbb_full_ht']:
                    if ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb','tthbb', 'ttbb', 'ttlf', 'ttcc','signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputH_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            #print(tmp)
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+['ttbb_full_ht'])
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)


        '''
        # Make vjets bkgs df
        dfList = []

        for i in vjets:
            for j in filein['columns'][f'{i}'].keys():
                tmp = []
                for var in pre_vars:
                    if ((var == 'process') & ('Jets' in j)):
                        inner_list = ['vjets'] * len(filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_ZH_M'].value.tolist())
                    elif ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
                tmp = np.transpose(np.asarray(tmp, dtype=object))
                print(tmp)
                tmpDF = pd.DataFrame(data=tmp, columns=pre_vars)
                dfList.append(tmpDF)
        vb_df = pd.concat(dfList, ignore_index=True)
        '''
        #s_df = s_df[s_df[genmatchreq] == True]
        #s_df = s_df.replace('old_ttZbb', 'ttZ')
        #s_df['sigVsbkg'] = s_df[['tthbb', 'ttzbb']].max(axis=1)
        #s_df['sigVsbkg'] = (s_df['ttzbb'] + s_df['tthbb'])/(s_df['ttzbb']+s_df['ttbb']+s_df['ttcc']+s_df['ttlf']+s_df['tthbb'])
        s_df['sigVsbkg'] = (s_df['ttzbb'] + s_df['tthbb'])/(s_df['ttzbb']+s_df['ttbb']+s_df['ttlf']+s_df['tthbb'])
        #s_df['sigVsbkg'] = s_df['signal']
        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]
        #b_df['sigVsbkg'] = (b_df['ttzbb'] + b_df['tthbb'])/(b_df['ttzbb']+b_df['ttbb']+b_df['ttcc']+b_df['ttlf']+b_df['tthbb'])
        #b_df['sigVsbkg'] = b_df[['tthbb', 'ttzbb']].max(axis=1)
        b_df['sigVsbkg'] = (b_df['ttzbb'] + b_df['tthbb'])/(b_df['ttzbb']+b_df['ttbb']+b_df['ttlf']+b_df['tthbb'])
        #b_df['sigVsbkg'] = b_df['signal']
        #b_df = pd.concat([b_df, vb_df])


        return s_df, b_df
        
    def get_NN_bins(self):
        cuts = dnn_cut(self.sb_df)
        genm = (self.sb_df['matchedGen_ZHbb_bb'] == True)
        # currently not binning by pt...
        nn_df = np.asarray(self.sb_df[ cuts & genm & (self.sb_df['ZH_pt'] < pthigh) & (self.sb_df['ZH_pt'] >= ptlow)]['sigVsbkg'].to_numpy(), dtype = float)
        norm_weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['ZH_pt'] < pthigh) & (self.sb_df['ZH_pt'] >= ptlow)]['norm_weight'].to_numpy(), dtype = float)
        weight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['ZH_pt'] < pthigh) & (self.sb_df['ZH_pt'] >= ptlow)]['genWeight'].to_numpy(), dtype = float)
        topptWeight = np.asarray(self.sb_df[cuts & genm & (self.sb_df['ZH_pt'] < pthigh) & (self.sb_df['ZH_pt'] >= ptlow)]['topptWeight'].to_numpy(), dtype = float)
        norm_weight = (topptWeight * norm_weight * np.sign(weight))
        quantiles = [0.0, .05, .25, .35, .50, .70, 1.0]

        nn_bins = weighted_quantile(nn_df, quantiles, norm_weight)
        nn_bins[0], nn_bins[-1] = 0,1
        #nn_bins[0] = 0
        print(nn_bins)
        #nn_bins = [0., 0.083, 0.431, 0.560, 0.736, 0.865, 1. ]
        #nn_bins = [0., 0.20, 0.40, 0.60, 0.80, 1.]
        #nn_bins = np.arange(0,1, 0.05)
        return nn_bins

    def plot_dnnHist(self, nn_bins):
        print(np.unique(self.sb_df.process, return_counts=True))
        cuts = dnn_cut(self.sb_df)
        sigs = ['ttZ', 'ttH']
        bkgs = ['tt_B', 'TTBar']#, 'vjets']


        label_dict = {
                "TTBar" : r'$\mathsf{t\bar{t}+\text{LF}}$, $\mathsf{t\bar{t}+c\bar{c}}$',
                "tt_B" : r'$\mathsf{t\bar{t}+b\bar{b}}$',
                "ttZ" : r'$\mathsf{t\bar{t}Z}$',
                "ttH" : r'$\mathsf{t\bar{t}H}$',
                "vjets": r'$\mathsf{V+jets}$'
        }
        
        for k in range(len(nn_bins) - 1, 0, -1):  # Iterate from the highest to the second highest bin
            lower_bound = nn_bins[k - 1]
            upper_bound = nn_bins[k]
            nn_cuts = (
                (self.sb_df['sigVsbkg'] >= lower_bound) &
                (self.sb_df['sigVsbkg'] < upper_bound)
            )
            #print(lower_bound)
            #print(upper_bound)
            nn_cuts = nn_cuts #& (self.sb_df['min_wpart_ZH_dR'] > 0.0)
            sumS = []
            sumB = []
            fig, ax = plt.subplots()
            fig.set_size_inches(3.75, 3.75)
            

            for i in bkgs+sigs:
                norm_weight = ((np.asarray(self.sb_df[cuts & nn_cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float))/(41.529))*137.596
                norm_weight2 = ((np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float))/(41.529))*137.596
                weight = np.asarray(self.sb_df[cuts & nn_cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
                weight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
                topptWeight = np.asarray(self.sb_df[cuts & nn_cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
                topptWeight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
    
                norm_weight = (topptWeight * norm_weight * np.sign(weight))
                norm_weight2 = (topptWeight2 * norm_weight2 * np.sign(weight2))
    
                if i in bkgs:
                    n, bins, patches = ax.hist(self.sb_df['ttbb_full_ht'][cuts & nn_cuts & (self.sb_df['process'] == i)], bins= 60,
                    #n, bins, patches = ax.hist(self.sb_df['ZH_M'][cuts & nn_cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5),
                        stacked=True, #bins= np.arange(0,0+5,1),
                        histtype='stepfilled', label=label_dict[i]) #weights=norm_weight)
                    #print("okay")
                #else:
                    #n, bins, patches = ax.hist(self.sb_df['ttbb_full_ht'][cuts & nn_cuts & (self.sb_df['process'] == i)], bins= 100,
                    #n, bins, patches = ax.hist(self.sb_df['ZH_M'][cuts & nn_cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5),
                        #stacked=False, #bins= np.arange(0,0+5,1),
                        #histtype='step', label=label_dict[i] + ' x10', weights=norm_weight*10)
                    #n, bins, patches = ax.hist(self.sb_df['ttbb_full_ht'][(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)], bins= 100,
                    #n, bins, patches = ax.hist(self.sb_df['ZH_M'][(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)], bins= np.arange(50,200+5,5),
                        #stacked=False, #bins= np.arange(0,0+5,1),
                        #histtype='step', label=label_dict[i] + '$_{GenMatch}$ x10', weights=norm_weight2*10)
                #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
                #        histtype='step', range= (0,1), label=f'{i}')
                ###################################################################################
                ##### FIXME going to add some functionality for yields and stat uncertainties #####
                ###################################################################################
                ht_cut = ((self.sb_df['ttbb_full_ht'] >= 750) & (self.sb_df['ttbb_full_ht'] < 7000))
                if i in sigs:
                    norm_weight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
                    weight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
                    topptWeight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & nn_cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
    
                    norm_weight_s_unc = (topptWeight_s_unc * norm_weight_s_unc * np.sign(weight_s_unc))
                    
                    s_yield = np.sum(norm_weight_s_unc)
                    s_uncertainty = np.sqrt(np.sum(norm_weight_s_unc**2))

                    print(i, "Matched yield:", s_yield, "+/-", s_uncertainty)

                    norm_weight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & nn_cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
                    weight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & nn_cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
                    topptWeight_s_unc = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & nn_cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
    
                    norm_weight_s_unc = (topptWeight_s_unc * norm_weight_s_unc * np.sign(weight_s_unc))
                    
                    s_yield = np.sum(norm_weight_s_unc)
                    s_uncertainty = np.sqrt(np.sum(norm_weight_s_unc**2))

                    print(i, "Unmatched yield:", s_yield, "+/-", s_uncertainty)
                if i in bkgs:
                    norm_weight_b_unc = np.asarray(self.sb_df[ht_cut & cuts & nn_cuts & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
                    weight_b_unc = np.asarray(self.sb_df[ht_cut & cuts & nn_cuts & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
                    topptWeight_b_unc = np.asarray(self.sb_df[ht_cut & cuts & nn_cuts & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
    
                    norm_weight_b_unc = (topptWeight_b_unc * norm_weight_b_unc * np.sign(weight_b_unc))
                    
                    b_yield = np.sum(norm_weight_b_unc)
                    #b_avg = (b_yield)/(len(norm_weight_b_unc))
                    #b_avg_unc = np.sqrt(b_avg**2 * (len(norm_weight_b_unc)))
                    #b4_avg = (b_yield)/(len(norm_weight_b_unc)*4)
                    #b4_avg_unc = np.sqrt(b4_avg**2 * (len(norm_weight_b_unc)*4))
                    b_uncertainty = np.sqrt(np.sum(norm_weight_b_unc**2))

                    print(i, "yield:", b_yield, "+/-", b_uncertainty)
                    print("n events:", len(norm_weight_b_unc))
                    #print(b_avg, b_avg_unc)
                    #print(b4_avg, b4_avg_unc)

                
                if i in sigs:
                    sumS.append(n)
                else:
                    sumB.append(n)
                #print(i, np.sum(n))
    
            #vals = pd.DataFrame(np.concatenate(sumS))
            #print(vals)
            #print(sum(np.sum(sumS,axis=0)), sum(np.sum(sumB,axis=0)))
            sumS = np.sum(sumS,axis=0)
            sumB = np.sum(sumB,axis=0)
            sob = sum(sumS)/sum(sumB)
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
            tex_x_corr = 0.42
            if ptrange == 'l':
                fig.text(tex_x_corr,0.59, r'200 < $\text{p}_T^\text{Z/H cand.}$ < 300 GeV', fontsize=10)
            if ptrange == 'm':
                fig.text(tex_x_corr,0.59, r'300 < $\text{p}_T^\text{Z/H cand.}$ < 450 GeV', fontsize=10)
            if ptrange == 'h':
                fig.text(tex_x_corr,0.59, r'$\text{p}_T^\text{Z/H cand.}$ > 450 GeV', fontsize=10)
            fig.text(tex_x_corr - .10,0.53, f'{lower_bound:4.3f} <= SPANet score < {upper_bound:4.4f}', fontsize=10)
            #fig.text(tex_x_corr,0.50, f'SoB: {sob:.2f}', fontsize=10)
    
            ax.set_ylabel("Events / 5 GeV", fontsize=10)
            #ax.set_ylabel("Counts", fontsize=10)
            #ax.set_xlabel(r"$m_{PNet}^\text{Z/H cand.}$", fontsize=10)
            ax.set_xlabel(r"$HT_\text{bb}$", fontsize=10)
            #ax.set_xlabel(r"$n~Topb~dR>1.2$", fontsize=10)
    
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=10)
            ax.tick_params(which='major', direction='in', top=True, right=True, labelsize=10, length=8)
            ax.tick_params(which='minor', direction='in', top=True, right=True, labelsize=10, length=4)
            #ax.set_yscale('log')
            #ax.set_xlim([bins[0],bins[-1]])
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=10, loc=0)
            fig.subplots_adjust(
                top=0.88,
                bottom=0.11,
                left=0.11,
                right=0.88,
                hspace=0.0,
                wspace=0.2)
            ax.set_ylim(0, 30)
    
            hep.cms.label("Work in progress", loc=0, ax=ax, fontsize=10)
            if ptrange == 'l':
                plt.savefig(f"pdf/sl_mass_{ptrange}_{dirname}_{k}.pdf", bbox_inches='tight')
            if ptrange == 'm':
                plt.savefig(f"pdf/sl_mass_{ptrange}_{dirname}_{k}.pdf", bbox_inches='tight')
            if ptrange == 'h':
                plt.savefig(f"pdf/sl_mass_{ptrange}_{dirname}_{k}.pdf", bbox_inches='tight')


if __name__ == '__main__':
    _ = DNN_datasets()
