import pandas as pd
import numpy as np

from coffea.util import load
import argparse
import outvars
import mplhep as hep

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

hep.style.use("CMS")
parser = argparse.ArgumentParser(description='Build datasets for NN training')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')


args=parser.parse_args()

filein = load(args.input)


def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] == 2) &
        (df_['ZH_bbvLscore'] >= 0.9105) &
        #(df_['ZH_bbvLscore'] >= 0.6) &
        (df_['n_ak4jets']   >= 5)              &
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
            (df_['ZH_pt'] > 200)
            #(df_['ZH_M'] > 75)       &
            #(df_['ZH_M'] < 145)      &
            #(df_['newgenm_NN'] <= 1) &
            #(df_['newgenm_NN'] > 0.)
    )
    return base_cuts



NN_vars = outvars.common_vars
sig_vars = outvars.common_vars+outvars.sig_vars
bkg_vars = outvars.common_vars+outvars.bkg_vars

#sig = ['ttHTobb__genMatch', 'ttHToNonbb__genMatch','TTZToBB__genMatch', 'TTZToQQ__genMatch', 'TTZToLLNuNu__genMatch']
#sig = ['ttHTobb__genMatch', 'TTZToBB__genMatch']
sig = ['ttHTobb__genMatch',
       #'ttHTobb__non_genMatch',
       'ttHToNonbb__genMatch',
       #'ttHToNonbb__non_genMatch',
       #'TTZToBB__genMatch',
       #'TTZToBB__non_genMatch',
       'TTZToQQ__genMatch',
       #'TTZToQQ__non_genMatch',
       #'TTZToLLNuNu__genMatch',
       #'TTZToLLNuNu__non_genMatch'
      ]
bkg = ["TTbb_Hadronic__tt+B",
       "TTbb_SemiLeptonic__tt+B",
       "TTbb_2L2Nu__tt+B",
       "TTToSemiLeptonic__tt+LF",
       "TTToSemiLeptonic__tt+C",
       "TTTo2L2Nu__tt+LF",
       "TTTo2L2Nu__tt+C",
       "TTToHadronic__tt+LF",
       "TTToHadronic__tt+C"]
       #"QCD_HT",
       #"WJetsToLNu_HT",
       #"DYJetsToLL_HT"]
       #]

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
    #dnn_vars = []
    #cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore', 'newgenm_NN','norm_weight', 'genWeight', 'topptWeight', 'tthbb']
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore','norm_weight', 'genWeight', 'topptWeight', 'signal']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df])
        self.nn_bins = self.get_NN_bins()
        self.plot_dnnHist()


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]

        genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')
        #print(genweight_df[0]['TTToSemiLeptonic__2017'])

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+[genmatchreq]:
                    if (var == 'norm_weight'):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'events_{var}'].value.tolist() # in latest Pocket Coffea add nominal
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb','tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'spanet_output_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+[genmatchreq])
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make bkg df
        '''
        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+['tt_type']:
                    if ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb','tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputH_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+['tt_type'])
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)
        '''
        for i in bkg:
            for j in filein['columns'][f'{i}'].keys():
                tmp = []
                for var in pre_vars:
                    if ((var == 'norm_weight') & ('TTBB' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    elif ((var == 'norm_weight') & ('TTBBtoLNu2Q' in j)):
                        inner_list = np.ones_like(filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'events_{var}'].value.tolist()) * 0.09
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'spanet_output_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
                tmp = np.transpose(np.asarray(tmp, dtype=object))
                #print(tmp)
                tmpDF = pd.DataFrame(data=tmp, columns=pre_vars)
                dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)
        b_df.loc[b_df['process'] == 'QCD_HT', 'process'] = 'QCD'
        b_df.loc[b_df['process'] == 'WJetsToLNu_HT', 'process'] = 'VJets'
        b_df.loc[b_df['process'] == 'DYJetsToLL_HT', 'process'] = 'VJets'

        #s_df = s_df[s_df[genmatchreq] == True]
        #s_df['sigVsbkg'] = (s_df['ttzbb'] + s_df['tthbb'])/(s_df['ttzbb']+s_df['ttbb']+s_df['ttlf']+s_df['tthbb'])
        s_df['sigVsbkg'] = s_df['signal']
        #s_df['sigVsbkg'] = (s_df['ttzbb'] + s_df['tthbb'])/(s_df['ttbb']+s_df['ttcc']+s_df['ttlf'])
        #s_df['sigVsbkg'] = s_df[['tthbb', 'ttzbb']].max(axis=1)
        #b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]
        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B') | (b_df['process'] == 'QCD') | (b_df['process'] == 'VJets')]
        #b_df['sigVsbkg'] = (b_df['ttzbb'] + b_df['tthbb'])/(b_df['ttzbb']+b_df['ttbb']+b_df['ttlf']+b_df['tthbb'])
        b_df['sigVsbkg'] = b_df['signal']
        b_df[genmatchreq] = False
        #b_df['sigVsbkg'] = (b_df['ttzbb'] + b_df['tthbb'])/(b_df['ttbb']+b_df['ttcc']+b_df['ttlf'])
        #b_df['sigVsbkg'] = b_df[['tthbb', 'ttzbb']].max(axis=1)
        return s_df, b_df

    def plot_dnnHist(self):
        fig, (ax, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        #fig.set_size_inches(3.75, 4.5)
        #fig, ax = plt.subplots(1,1)
        fig.set_size_inches(3.75, 3.75)
        print(np.unique(self.sb_df.process, return_counts=True))
        cuts = dnn_cut(self.sb_df)
        plotcut = plot_cut(self.sb_df)
        sig = ['old_ttZbb', 'ttH']
        bkg = ['tt_B', 'TTBar', 'VJets']
        sumS = []
        sumB = []

        label_dict = {
                "TTBar" : r'$\mathsf{t\bar{t}+\text{LF}}$, $\mathsf{t\bar{t}+c\bar{c}}$',
                "tt_B" : r'$\mathsf{t\bar{t}+b\bar{b}}$',
                "old_ttZbb" : r'$\mathsf{t\bar{t}Z}$',
                "ttH" : r'$\mathsf{t\bar{t}H}$',
                "VJets": r'$\mathsf{V+jets}$',
                "QCD": r'$\mathsf{QCD}$'
        }

        for i in sig + bkg:
            norm_weight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
            print(i, np.unique(self.sb_df[(self.sb_df['process'] == i)]['norm_weight'],return_counts=True))
            weight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            topptWeight = np.asarray(self.sb_df[cuts & plotcut & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)

            norm_weight = (topptWeight * norm_weight * np.sign(weight))

            norm_weight = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
            norm_weight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & plotcut & (self.sb_df['process'] == i)]['norm_weight'].to_numpy(), dtype = float)
            weight = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            weight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & plotcut & (self.sb_df['process'] == i)]['genWeight'].to_numpy(), dtype = float)
            topptWeight = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
            topptWeight2 = np.asarray(self.sb_df[(self.sb_df[genmatchreq] == True) & cuts & plotcut & (self.sb_df['process'] == i)]['topptWeight'].to_numpy(), dtype = float)
            norm_weight = (topptWeight * norm_weight * np.sign(weight))
            norm_weight2 = (topptWeight2 * norm_weight2 * np.sign(weight2))

            if i in bkg:
                print(i, norm_weight)
                n, bins, patches = ax.hist(self.sb_df['sigVsbkg'][(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
                    histtype='step', range= (0,1), label=label_dict[i], weights=norm_weight, density=False)
            else:
                n, bins, patches = ax.hist(self.sb_df['sigVsbkg'][(self.sb_df[genmatchreq] == True) & cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
            #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=10, stacked=False,
                        histtype='step', range= (0,1), label=label_dict[i] + "$_{matched}$", weights=norm_weight2, density=False)
                #n, bins, patches = ax.hist(self.sb_df['sigVsbkg'][(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
            #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=10, stacked=False,
                        #histtype='step', range= (0,1), label=label_dict[i], weights=norm_weight, density=True)

            #n, bins, patches = ax.hist(self.sb_df['sigVsbkg'][cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
            #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=10, stacked=False,
                  #  histtype='step', range= (0,1), label=label_dict[i], weights=norm_weight)

            #n, bins, patches = ax.hist(self.sb_df['newgenm_NN'][cuts & plotcut & (self.sb_df['process'] == i)], bins=self.nn_bins, stacked=False,
            #        histtype='step', range= (0,1), label=f'{i}')
            if i in sig:
                sumS.append(n)
            elif i in bkg:
                sumB.append(n)
            print(i, np.sum(n))
            print("N QCD", self.sb_df['sigVsbkg'][(self.sb_df[genmatchreq] == False) & cuts & plotcut & (self.sb_df['process'] == "QCD")])

        #vals = pd.DataFrame(np.concatenate(sumS))
        #print(vals)
        print(np.sum(sumS,axis=0), np.sum(sumB,axis=0))
        sumS = np.sum(sumS,axis=0)
        sumB = np.sum(sumB,axis=0)
        bin_c = (bins[1:]+bins[:-1])/2
        print("SoSqrtB", sumS/np.sqrt(sumB))
        ax2.errorbar(x=bin_c, y = sumS/np.sqrt(sumB), xerr=(bins[1:]-bins[:-1])/2,
                fmt='.', color='k', label=r'S/$\sqrt{\mathrm{B}}$')
        ax2.tick_params(which='both', direction='in', top=True, right=True, labelsize=10)
        ax2.set_ylabel(r'$\mathrm{S/}\sqrt{\mathrm{B}}$', fontsize=10)
        #ax2.yaxis.set_label_coords(-0.12,0.45)
        ax.set_xlabel(r'SPANet Signal Score', fontsize=10)
        #ax.set_ylabel(r'Density', fontsize=10)
        ax2.grid(True)
        ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=10)
        ax.set_yscale('log')
        ax.set_xlim([bins[0],bins[-1]])
        ax.set_xlim([0.8,bins[-1]])
        #ax.set_ylim(10**-2,30)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=10)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0,
            wspace=0.2)
        hep.cms.label(exp = "", label="Private work (CMS simulation)", loc=0, ax=ax, fontsize=10, data=True, rlabel='(13.6 TeV)')
        #hep.cms.lumitext("41.48")

        plt.savefig("sl_SPANET_2024_SoB.pdf", bbox_inches='tight')

    def get_NN_bins(self):
        cuts = dnn_cut(self.sb_df)
        genm = (self.sb_df[genmatchreq] == True)
        # currently not binning by pt...
        nn_df = np.asarray(self.sb_df[ cuts & genm]['sigVsbkg'].to_numpy(), dtype = float)
        norm_weight = np.asarray(self.sb_df[cuts & genm]['norm_weight'].to_numpy(), dtype = float)
        weight = np.asarray(self.sb_df[cuts & genm]['genWeight'].to_numpy(), dtype = float)
        topptWeight = np.asarray(self.sb_df[cuts & genm]['topptWeight'].to_numpy(), dtype = float)
        norm_weight = (topptWeight * norm_weight * np.sign(weight))
        quantiles = [0.0, .05, .25, .35, .50, .70, 1.0]
        #quantiles = [0.0, .25, .5, .75, 1.0]

        nn_bins = weighted_quantile(nn_df, quantiles, sample_weight=norm_weight)
        #nn_bins = weighted_quantile(nn_df, quantiles)
        nn_bins[0], nn_bins[-1] = 0,1
        #nn_bins[0] = 0
        print(nn_bins)
        #nn_bins = [0., 0.083, 0.431, 0.560, 0.736, 0.865, 1. ]
        #nn_bins = [0., 0.20, 0.40, 0.60, 0.80, 1.]
        #nn_bins = np.arange(0,1.05, 0.05)
        return nn_bins




if __name__ == '__main__':
    _ = DNN_datasets()
