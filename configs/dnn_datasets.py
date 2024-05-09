import pandas as pd
import numpy as np
import os
import sys
import awkward as ak
from glob import glob
from keras.utils import to_categorical
from sklearn.preprocessing import LabelEncoder

from coffea.util import load
import argparse
import outvars

parser = argparse.ArgumentParser(description='Build datasets for NN training')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] >= 2) &
        (df_['ZH_bbvLscore'] >= 0.8) &
        (df_['n_ak4jets']   >= 5)             &
        #( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        #(df_['passNotHadLep'] == 1) & # might add
        (df_['ZH_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['ZH_M']        >= 50)            &
        (df_['ZH_M']        <= 200)
    )
    return base_cuts

NN_vars = outvars.NN_vars
sig_vars = outvars.NN_vars+outvars.sig_vars
bkg_vars = outvars.NN_vars+outvars.bkg_vars

sig = ['ttHTobb__genMatch', 'ttHToNonbb__genMatch','TTZToBB__genMatch', 'TTZToQQ__genMatch', 'TTZToLLNuNu__genMatch']
#sig = ['ttHTobb', 'TTZToBB']
bkg = [
    "TTbb_SemiLeptonic__TTbbSemiLeptonic_tt+B",
    "TTbb_SemiLeptonic__TTbbSemiLeptonic_tt+LF",
    "TTbb_SemiLeptonic__TTbbSemiLeptonic_tt+C",
    "TTToSemiLeptonic__TTToSemiLeptonic_tt+LF",
    "TTToSemiLeptonic__TTToSemiLeptonic_tt+C",
    "TTToSemiLeptonic__TTToSemiLeptonic_tt+B"
    ]
genmatchreq = 'matchedGen_ZHbb_bb'

class DNN_datasets:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    #bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = ['process','ZH_pt','MET_pt','ZH_M', 'ZH_bbvLscore']
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        #print(self.s_df[self.s_df['matchedGen_ZHbb_bb'] == True])
        self.prep_class()
        self.sep_test_train()


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]

        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+[genmatchreq]:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+[genmatchreq])
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make signal df

        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in pre_vars+['tt_type']:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=pre_vars+['tt_type'])
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        s_df = s_df[s_df[genmatchreq] == True]
        b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

        return s_df, b_df

    def prep_class(self):
        self.s_df['label'] = 2
        self.b_df['label'] = np.where(self.b_df['process'] == 'TTBar', 0, 1)
        del self.s_df[genmatchreq], self.b_df['tt_type']
        #del self.b_df[genmatchreq], self.s_df['tt_type']
        print(self.s_df)
        sb_df = pd.concat([self.s_df,self.b_df])
        #
        #sb_df.drop(columns="process", inplace=True)
        #sb_df = sb_df.astype(np.float64)
        #print(sb_df.columns.to_series()[np.isinf(sb_df).any()])
        #print(sb_df.index[np.isinf(sb_df).any(1)])
        #print(sb_df.loc[139, sb_df.columns.to_series()[np.isinf(sb_df).any()]])
        #
        sb_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        sb_df.dropna(how="any", inplace=True)
        print(sb_df)
        print(sum(sb_df['process'] == 'sig'))
        encoder = LabelEncoder()
        encoder.fit(sb_df['label'])
        encoded_labels = encoder.transform(sb_df['label'])
        onehot_labels = to_categorical(encoded_labels)
        sb_df['label'] = onehot_labels.tolist()
        self.sb_df = sb_df[dnn_cut(sb_df)].sample(frac=1).reset_index(drop=True) # to shuffle dataset
        print(np.unique(self.sb_df['label'],return_counts=True))
        for v in self.cut_vars:
            if v not in self.dnn_vars:
                del self.sb_df[v]

    def sep_test_train(self):
        train_df = self.sb_df.sample(frac=.75, random_state=1)
        test_df  = self.sb_df.drop(train_df.index).copy().sample(frac=1).reset_index(drop=True)
        train_df.to_pickle('./trainXY.pkl')
        test_df.to_pickle('./testXY.pkl')

if __name__ == '__main__':
    _ = DNN_datasets()
