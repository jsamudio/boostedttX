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

sig = ['ttHTobb__genMatch',
       'ttHToNonbb__genMatch',
       'TTZToBB__genMatch',
       'TTZToQQ__genMatch',
       'TTZToLLNuNu__genMatch']

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
processes = ['ttZ', 'ttH', 'TTBar', 'ttbb']


class outputToPkl:

    def __init__(self):

        self.filein = filein
        self.sig_vars = outvars.NN_vars+outvars.sig_vars+outvars.weight_vars
        self.bkg_vars = outvars.NN_vars+outvars.bkg_vars+outvars.weight_vars

        self.s_df, self.b_df = self.get_sigbkg()
        self.exportPkl()

        print(self.s_df)
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
            for j in self.filein['columns'][f'{i}'].keys():
                for var in self.sig_vars:
                    if (var == 'norm_weight'):
                        inner_list = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    else:
                        inner_list = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=self.sig_vars)
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make bkg df

        for i in bkg:
            tmp = []
            for j in self.filein['columns'][f'{i}'].keys():
                for var in self.bkg_vars:
                    if ((var == 'norm_weight') & ('TTbb' not in j)):
                        inner_list = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                        inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                    else:
                        inner_list = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=self.bkg_vars)
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

        '''
        Background is automatically trimmed to just TTBar and ttbb, and if genmatched
        signal is needed this can be done after this function.
        '''

        return s_df, b_df
    
    def sepByProcess(self):
    
    def exportPkl(self):
        self.s_df.to_pickle("./dummySig.pkl")

if __name__ == '__main__':
    _ = outputToPkl()