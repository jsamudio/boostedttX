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
       #'TTZToBB__genMatch',
       'TTZToQQ__genMatch',
       #'TTZToLLNuNu__genMatch',
       'ttHTobb__non_genMatch',
       'ttHToNonbb__non_genMatch',
       #'TTZToBB__non_genMatch',
       'TTZToQQ__non_genMatch']
       #'TTZToLLNuNu__non_genMatch']

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


# How to keep these centralized and permanently updated?
processes = ['ttZ', 'ttH', 'TTBar', 'tt_B', 'QCD', 'VJets']


class outputToPkl:

    def __init__(self):

        self.filein = filein
        #print(filein['columns'])
        #self.sig_vars = outvars.NN_vars+outvars.sig_vars+outvars.weight_vars+['newgenm_NN'] #['ttzbb', 'tthbb', 'genZHpt']
        self.sig_vars = outvars.common_vars+outvars.sig_vars+outvars.weight_vars+['signal']
        #self.bkg_vars = outvars.common_vars+outvars.bkg_vars+outvars.weight_vars+['signal']
        self.bkg_vars = outvars.common_vars+outvars.weight_vars+['signal']

        self.s_df, self.b_df, self.data_obs = self.get_sigbkg()
        self.sb_df = pd.concat([self.s_df,self.b_df,self.data_obs])
        self.sb_df.loc[self.sb_df['process'] == 'old_ttZbb', 'process'] = 'ttZ'
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
        print('-' * 60)
        print('WRITING TO PKL (OPTIMIZED)')
        print('-' * 60)
        
        genweight_df = pd.DataFrame.from_dict(self.filein['sum_signOf_genweights'], orient='index')

        # ---------------------------------------------------------
        # 1. Build Signal DataFrame
        # ---------------------------------------------------------
        dfList = []
        for i in sig:
            for j in self.filein['columns'][f'{i}'].keys():
                # Fetch genweight and base dictionary ONCE per dataset
                genweight = genweight_df[0][f'{j}']
                base_dict = self.filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal']
                
                tmp_data = {}
                for var in self.sig_vars:
                    if var == 'norm_weight':
                        # Use Numpy array for fast, vectorized division
                        arr = np.array(base_dict[f'events_{var}'].value)
                        tmp_data[var] = arr / genweight
                    elif var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']:
                        tmp_data[var] = base_dict[f'spanet_output_{var}'].value
                    else:
                        tmp_data[var] = base_dict[f'events_{var}'].value
                
                # Create DataFrame directly from the dictionary (avoids transpose)
                dfList.append(pd.DataFrame(tmp_data))
                
        s_df = pd.concat(dfList, ignore_index=True)

        # ---------------------------------------------------------
        # 2. Build Background DataFrame
        # ---------------------------------------------------------
        ttbb_xsecs = {
            "TTbb_SemiLeptonic__tt+B": 17.36875,
            "TTbb_2L2Nu__tt+B": 3.74766,
            "TTbb_Hadronic__tt+B": 19.15297
        }
        lumi_pb = 109.95 * 1000 # Convert fb-1 to pb-1

        dfList = []
        for i in bkg:
            for j in self.filein['columns'][f'{i}'].keys():
                genweight = genweight_df[0][f'{j}'] # This is your sum_of_genweights
                base_dict = self.filein['columns'][f'{i}'][f'{j}']['btag_mask']['nominal']
                
                tmp_data = {}
                for var in self.bkg_vars:
                    if var == 'norm_weight':
                        arr = np.array(base_dict[f'events_{var}'].value)
                        
                        # Apply new xsec * lumi formula specifically for TTbb samples
                        if i in ttbb_xsecs:
                            # Formula: (event_weight / sum_genweights) * xsec * lumi
                            tmp_data[var] = (arr / genweight) * (ttbb_xsecs[i] * lumi_pb)
                        else:
                            # Standard normalization for all other backgrounds
                            tmp_data[var] = arr / genweight
                            
                    elif var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']:
                        tmp_data[var] = base_dict[f'spanet_output_{var}'].value
                    else:
                        tmp_data[var] = base_dict[f'events_{var}'].value
                        
                dfList.append(pd.DataFrame(tmp_data))
                
        b_df = pd.concat(dfList, ignore_index=True)

        # ---------------------------------------------------------
        # 3. Clean up and rename background processes
        # ---------------------------------------------------------
        # Use .replace() for much faster bulk renaming
        rename_map = {
            'QCD_HT': 'QCD',
            'WJetsToLNu_HT': 'VJets',
            'DYJetsToLL_HT': 'VJets'
        }
        b_df['process'] = b_df['process'].replace(rename_map)

        # Filter background to strictly required processes
        allowed_bkgs = ['TTBar', 'tt_B', 'QCD', 'VJets']
        b_df = b_df[b_df['process'].isin(allowed_bkgs)]

        # ---------------------------------------------------------
        # 4. Mock Data Obs
        # ---------------------------------------------------------
        data_df = s_df.copy() # FIXME once we have data it will go here
        data_df = data_df.assign(process='data_obs')

        print('-' * 60)
        print('Done')
        print('-' * 60)

        return s_df, b_df, data_df
    
    def exportByProcess(self):
        for i in processes + ['data_obs']:
            export_df = self.sb_df[(self.sb_df['process'] == i)]
            export_df.to_pickle(f"pickled/{dirname}_{i}.pkl")
            

if __name__ == '__main__':
    _ = outputToPkl()