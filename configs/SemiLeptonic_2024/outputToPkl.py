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

# Master lists of allowed samples.
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
       "TTToHadronic__tt+C",
       #"QCD_HT",
       "WJets",
       "DYJets"]

processes = ['ttZ', 'ttH', 'TTBar', 'tt_B', 'QCD', 'VJets']

class outputToPkl:

    def __init__(self):
        self.filein = filein
        
        # Use dict.fromkeys to safely remove duplicate column names while preserving order
        self.sig_vars = list(dict.fromkeys(outvars.common_vars + outvars.sig_vars + outvars.weight_vars + ['signal']))
        self.bkg_vars = list(dict.fromkeys(outvars.common_vars + outvars.weight_vars + ['signal']))

        # 1. Dynamically discover all variations in this specific .coffea file
        available_sigs = [s for s in sig if s in self.filein['columns']]
        available_bkgs = [b for b in bkg if b in self.filein['columns']]
        
        if available_sigs:
            first_sample = available_sigs[0]
        elif available_bkgs:
            first_sample = available_bkgs[0]
        else:
            raise KeyError("Neither allowed signal nor background samples were found in the input file.")

        first_dataset = list(self.filein['columns'][first_sample].keys())[0]
        available_variations = list(self.filein['columns'][first_sample][first_dataset]['btag_mask'].keys())
        
        print(f"Found variations in file: {available_variations}")

        # 2. Loop through every variation and create a dedicated .pkl
        for var in available_variations:
            self.s_df, self.b_df, self.data_obs = self.get_sigbkg(variation=var)
            self.sb_df = pd.concat([self.s_df, self.b_df, self.data_obs])
            self.sb_df.loc[self.sb_df['process'] == 'old_ttZbb', 'process'] = 'ttZ'
            self.sb_df['nnscore'] = self.sb_df['signal']
            
            self.exportByProcess(variation=var)
        
    def get_sigbkg(self, variation='nominal'):
        print('-' * 60)
        print(f'WRITING TO PKL (OPTIMIZED) - VARIATION: {variation}')
        print('-' * 60)
        
        genweight_df = pd.DataFrame.from_dict(self.filein['sum_signOf_genweights'], orient='index')

        # --- Signal ---
        dfList = []
        for i in sig:
            if i not in self.filein['columns']:
                print(f"Skipping {i} (not found in input file)")
                continue
                
            for j in self.filein['columns'][f'{i}'].keys():
                genweight = genweight_df[0][f'{j}']
                base_dict = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][variation]
                
                tmp_data = {}
                for var in self.sig_vars:
                    if var == 'norm_weight':
                        arr = np.array(base_dict[f'events_{var}'].value)
                        tmp_data[var] = arr / genweight
                    elif var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']:
                        tmp_data[var] = base_dict[f'spanet_output_{var}'].value
                    else:
                        tmp_data[var] = base_dict[f'events_{var}'].value
                
                dfList.append(pd.DataFrame(tmp_data))
                
        if dfList:
            s_df = pd.concat(dfList, ignore_index=True)
        else:
            s_df = pd.DataFrame(columns=self.sig_vars)

        dfList = []
        for i in bkg:
            if i not in self.filein['columns']:
                print(f"Skipping {i} (not found in input file)")
                continue
                
            for j in self.filein['columns'][f'{i}'].keys():
                genweight = genweight_df[0][f'{j}'] 
                base_dict = self.filein['columns'][f'{i}'][f'{j}']['btag_mask'][variation]
                
                tmp_data = {}
                for var in self.bkg_vars:
                    if var == 'norm_weight':
                        arr = np.array(base_dict[f'events_{var}'].value)
                        tmp_data[var] = arr / genweight
                    elif var in ['ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc', 'signal']:
                        tmp_data[var] = base_dict[f'spanet_output_{var}'].value
                    else:
                        tmp_data[var] = base_dict[f'events_{var}'].value
                        
                dfList.append(pd.DataFrame(tmp_data))
                
        if dfList:
            b_df = pd.concat(dfList, ignore_index=True)
        else:
            # ONLY use bkg_vars to avoid duplicate column generation
            b_df = pd.DataFrame(columns=self.bkg_vars)

        if not b_df.empty and 'process' in b_df.columns:
            rename_map = {
                'QCD_HT': 'QCD',
                'WJets': 'VJets',
                'DYJets': 'VJets'
            }
            b_df['process'] = b_df['process'].replace(rename_map)
    
            b_df.loc[b_df['process'].str.contains(r'tt\+B', na=False), 'process'] = 'tt_B'
            b_df.loc[b_df['process'].str.contains(r'tt\+LF|tt\+C', na=False), 'process'] = 'TTBar'
    
            allowed_bkgs = ['TTBar', 'tt_B', 'QCD', 'VJets']
            b_df = b_df[b_df['process'].isin(allowed_bkgs)]

        # Mock Data so no accidents happen
        data_df = s_df.copy() 
        if not data_df.empty:
            data_df = data_df.assign(process='data_obs')

        return s_df, b_df, data_df
    
    def exportByProcess(self, variation):
        formatted_var = variation.replace("_up", "Up").replace("_down", "Down")
        var_suffix = f"_{formatted_var}" if formatted_var != "nominal" else ""
        
        for i in processes + ['data_obs']:
            if i == 'data_obs' and variation != 'nominal': 
                continue
                
            export_df = self.sb_df[(self.sb_df['process'] == i)]
            
            if export_df.empty:
                continue
                
            output_name = f"pickled/Inference_{i}{var_suffix}.pkl"
            os.makedirs(os.path.dirname(output_name), exist_ok=True)
            
            export_df.to_pickle(output_name)
            print(f"Saved {output_name}")
            

if __name__ == '__main__':
    _ = outputToPkl()