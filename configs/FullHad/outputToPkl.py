import pandas as pd
import numpy as np
import os
import sys
from glob import glob

from coffea.util import load
import argparse
import outvars

parser = argparse.ArgumentParser(description='Change .coffea output to .pkl')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)
dirname = os.path.dirname(args.input)

sig = ['ttHTobb__genMatch',
       'ttHToNonbb__genMatch',
       'TTZToBB__genMatch',
       'TTZToQQ__genMatch',
       'TTZToLLNuNu__genMatch',
       'ttHTobb__non_genMatch',
       'ttHToNonbb__non_genMatch',
       'TTZToBB__non_genMatch',
       'TTZToQQ__non_genMatch',
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


# How to keep these centralized and permanently updated?
processes = ['ttZ', 'ttH', 'TTBar', 'tt_B']


class outputToPkl:

    def __init__(self):

        self.filein = filein
        #print(filein['columns'])
        #self.sig_vars = outvars.NN_vars+outvars.sig_vars+outvars.weight_vars+['newgenm_NN'] #['ttzbb', 'tthbb', 'genZHpt']
        self.sig_vars = outvars.NN_vars+outvars.sig_vars+outvars.weight_vars+['signal']
        self.bkg_vars = outvars.NN_vars+outvars.bkg_vars+outvars.weight_vars+['signal']

        self.s_df, self.b_df, self.data_obs = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df,self.data_obs])
        #self.sb_df['nnscore'] = (self.sb_df['ttzbb'] + self.sb_df['tthbb'])/(self.sb_df['ttzbb']+self.sb_df['ttbb']+self.sb_df['ttcc']+self.sb_df['ttlf']+self.sb_df['tthbb'])
        # FIXME this is for the nottcc only
        self.sb_df['nnscore'] = self.sb_df['signal']
        #self.sb_df['nnscore'] = (self.sb_df['ttzbb'] + self.sb_df['tthbb'])/(self.sb_df['ttzbb']+self.sb_df['ttbb']+self.sb_df['ttlf']+self.sb_df['tthbb'])
        #self.sb_df['nnscore'] = self.sb_df[['tthbb', 'ttzbb']].max(axis=1)
        self.exportByProcess()

        '''
        psuedo code like:

        get_sigbkg(self)
        sepProcess(self) # this would break down the df by process (tagged MC process specifically ttZ, ttH, TTBar, ttW, etc.)
        exportPkl(self) # this should export and save a pkl for each tagged MC process 
        '''
        
    def get_sigbkg(self):
        genweight_df = pd.DataFrame.from_dict(self.filein['sum_signOf_genweights'], orient='index')

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            print(filein['columns'].keys())
            for j in filein['columns'][f'{i}'].keys():
                for var in self.sig_vars:
                    if (var == 'norm_weight'):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputH_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=self.sig_vars)
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make bkg df

        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in self.bkg_vars:
                    if ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    #elif (var in ['ttzbb']):
                    #    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputZ_{var}'].value.tolist()
                    elif (var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']):
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'spanet_outputH_{var}'].value.tolist()
                    else:
                        inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            #print(tmp)
            tmpDF = pd.DataFrame(data=tmp, columns=self.bkg_vars)
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

        data_df = s_df.copy()
        data_df = data_df.assign(process = 'data_obs')
        #print(s_df)

        '''
        Background is automatically trimmed to just TTBar and ttbb, and if genmatched
        signal is needed this can be done after this function.
        '''

        return s_df, b_df, data_df
    
    def exportByProcess(self):
        for i in processes + ['data_obs']:
            export_df = self.sb_df[(self.sb_df['process'] == i)]
            export_df.to_pickle(f"pickled/{dirname}_{i}.pkl")
            

if __name__ == '__main__':
    _ = outputToPkl()