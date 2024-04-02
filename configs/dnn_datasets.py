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

parser = argparse.ArgumentParser(description='Build datasets for NN training')

parser.add_argument('--input', '-i', type=str, help='Input .coffea file')

args=parser.parse_args()

filein = load(args.input)

def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] >= 2) &
        (df_['ZH_bbvLscore'] >= 0.8)
        )
    return base_cuts

NN_vars = [
    'matchedGen_ZHbb_bb', 'process', 'tt_type',
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score','outZH_b2_score',
    'outZH_q1_pt','outZH_q2_pt',
    'outZH_q1_score','outZH_q2_score',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1','outZH_q_q_dr_nearb2',
    'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    'outZH_b1_qq_dr','outZH_b2_qq_dr',
    'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr','ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    'l_b2_mtb',
    #
    'ZH_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_ZHbb',
    'outZH_max_ak8pnetMass',
    'outZH_b12_m', 'outZH_b12_dr',
    'ht_b', 'ht_outZH',
    #
    'ak4_bestb_inZH',
    'ak4_worstb_inZH',
    #
    'nonZHbb_q1_dr',
    'nonZHbb_b1_dr',
    'inZHb_outZHb_dr',
    #
    'ZH_l_dr', 'ZH_l_invM',
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZH', 'n_q_inZH',
    'n_b_outZH', 'n_q_outZH', "ZH_bbvLscore"]

#sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']
sig = ['ttHTobb', 'TTZToBB']
bkg = ["TTbb_SemiLeptonic", "TTToSemiLeptonic"]
genmatchreq = 'matchedGen_ZHbb_bb'

class DNN_datasets:
    #Prepare datasets for MVA training

    #sig = ['ttH', 'ttZ']
    bkg = ['TTBar', 'ttbb']
    dnn_vars = NN_vars
    cut_vars = []
    output_dir = './nn_files'

    def __init__(self):
        self.s_df, self.b_df = self.get_sigbkg()
        print(self.s_df[self.s_df['matchedGen_ZHbb_bb'] == True])
        self.prep_class()
        self.sep_test_train()


    def get_sigbkg(self):
        pre_vars = self.dnn_vars + [v for v in self.cut_vars if v not in self.dnn_vars]



        dfList = []

        #Make signal df

        for i in sig:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in NN_vars:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=NN_vars)
            dfList.append(tmpDF)
        s_df = pd.concat(dfList, ignore_index=True)

        dfList = []

        #Make signal df

        for i in bkg:
            tmp = []
            for j in filein['columns'][f'{i}'].keys():
                for var in NN_vars:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    tmp.append(inner_list)
            tmp = np.transpose(np.asarray(tmp, dtype=object))
            tmpDF = pd.DataFrame(data=tmp, columns=NN_vars)
            dfList.append(tmpDF)
        b_df = pd.concat(dfList, ignore_index=True)

        s_df = s_df[s_df[genmatchreq] == True]

        return s_df, b_df

    def prep_class(self):
        self.s_df['label'] = 2
        self.b_df['label'] = np.where(self.b_df['process'] == 'TTbar', 0, 1)
        del self.s_df[genmatchreq], self.b_df['tt_type']
        del self.b_df[genmatchreq], self.s_df['tt_type']
        sb_df = pd.concat([self.s_df,self.b_df])
        encoder = LabelEncoder()
        encoder.fit(sb_df['label'])
        encoded_labels = encoder.transform(sb_df['label'])
        onehot_labels = to_categorical(encoded_labels)
        sb_df['label'] = onehot_labels.tolist()
        self.sb_df = sb_df[dnn_cut(sb_df)].sample(frac=1).reset_index(drop=True) # to shuffle dataset
        self.sb_df.fillna(value=0, inplace=True)
        print(np.unique(self.sb_df['process'],return_counts=True))
        del self.sb_df['process']
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
