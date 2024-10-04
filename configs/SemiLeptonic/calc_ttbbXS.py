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

lumi = 41.48

class calc_ttbbXS:

    def __init__(self):
        self.s_df, self.b_df = self.get_ttBar()
        self.calc_ttbbXS()

    def get_ttBar(self):
        #ttBar = ['TTToHadronic__tt+B', 'TTToSemiLeptonic__tt+B', 'TTTo2L2Nu__tt+B']
        #ttBar = ['TTToSemiLeptonic__tt+B']
        #ttbb = ['TTbb_SemiLeptonic__tt+B']
        #ttBar = ['TTToHadronic__tt+B']
        #ttbb = ['TTbb_Hadronic__tt+B']
        ttBar = ['TTTo2L2Nu__tt+B']
        ttbb = ['TTbb_2L2Nu__tt+B']
        pre_vars = ['genWeight']

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

        for i in ttbb:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars)
            dfList.append(tmpDF)

        b_df = pd.concat(dfList, ignore_index=True)

        return s_df, b_df

    def calc_ttbbXS(self):
        # pull in xsec from metadata? or just hard code it
        # next, next get sums:
        #xs = 364.018
        #xs = 380.095
        xs = 88.29

        #sum_signOf_genW = filein['sum_signOf_genweights']['TTToSemiLeptonic__2017']
        #sum_signOf_genW = filein['sum_signOf_genweights']['TTToHadronic__2017']
        sum_signOf_genW = filein['sum_signOf_genweights']['TTTo2L2Nu__2017']
        sumW = sum(np.sign(self.s_df['genWeight']))
        ttbbXS = xs * (sumW/sum_signOf_genW)
        print("The yield XS from 5FS:", ttbbXS)

        ttbbSumW = sum(np.sign(self.b_df['genWeight']))
        #ttbbSum_signOf_genW = filein['sum_signOf_genweights']['TTbb_SemiLeptonic__2017']
        #ttbbSum_signOf_genW = filein['sum_signOf_genweights']['TTbb_Hadronic__2017']
        ttbbSum_signOf_genW = filein['sum_signOf_genweights']['TTbb_2L2Nu__2017']
        invertRatio = ttbbSum_signOf_genW/ttbbSumW
        print("Total xsec for 4FS:", ttbbXS*invertRatio)



if __name__ == '__main__':
    _ = calc_ttbbXS()

