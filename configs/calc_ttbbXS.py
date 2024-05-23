import pandas as pd
import numpy as np
import os
import sys
import awkward as ak
from coffea.util import load

import argparse

parser = argparse.ArgumentParser(description='Calculate 4FS xsec from 5FS')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def get_ttBar(self):
    ttBar = ['TTToHadronic__tt+B', 'TTToSemiLeptonic__tt+B', 'TTTo2L2Nu__tt+B']
    pre_vars = ['tt_B', 'lumi', 'genWeight', 'sum_sign_genw']

    dfList = []

    #Make input df

    for i in ttBar:
        tmp = []
        for j in filein['columns'][f'{i}'].keys():
            for var in pre_vars:
                inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                tmp.append(inner_list)
        tmp = np.transpose(np.asarray(tmp, dtype=object))
        tmpDF = pd.DataFrame(data=tmp, columns=pre_vars)
        dfList.append(tmpDF)
    s_df = pd.concat(dfList, ignore_index=True)

    dfList = []

    return s_df

def calc_ttbbXS(df):
    # pull in xsec from metadata? or just hard code it
    # next, next get sums:

    #FIXME
    sumW = s_df['genWeight'][(s_df['tt_B'] == True) & (isHadronic??)]
