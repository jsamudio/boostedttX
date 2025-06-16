#                    #
##                  ##
######################
### Build datacard ###
######################
######################

import numpy as np
import pandas as pd
import json
import os
import sys
import functools

def getZhbbWeight(df_, year):
    tot_weight = (((df_['norm_weight']/(41.529))*(137.596+62+120)) * np.sign(df_['genWeight']) * df_['topptWeight']) # and other weights
    return tot_weight

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
    values = np.array(values, dtype='float')
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight, dtype='float')
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

class DataCardShapes():
    '''
    Handle format of datacard shapes
    '''
    years = ['2017'] # more later, centralized
    file_dir = './pickled'
    ref_samples = ['ttZ', 'ttH']
    nn = 'nnscore'
    hist_dict = {}
    
    def __init__(self, recopt_bins, recosdM_bins, n_NN_bins=10, nn=nn, isblind=True): # will have to change isblind for final fit
        self.pt_bins = recopt_bins+[500]
        self.sdM_bins = recosdM_bins
        self.nn = nn
        self.isblind = isblind
        self.nn_bins = None # normally a dict by year
        self.init_hist_funcs()
        
    def init_hist_funcs(self):
        for y in  self.years:
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/SpanetInferenceTruncated_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/SpanetInferenceDoubleWithGenPt_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/SpanetInferenceAssignment_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/SpanetBalanced24SDM_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/SpanetNoTTCC_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/DNNInference_{s}.pkl'))
            #get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/XbbVsQCD_{s}.pkl'))
            get_pickle= (lambda s: pd.read_pickle(f'{self.file_dir}/USCMSposter_inference_{s}.pkl'))
            df = pd.concat([get_pickle(s) for s in self.ref_samples], axis='rows', ignore_index=True)
            #print(df)
            #
            #print(self.nn)
            df = df[((df[self.nn]>=0.0) & (df['matchedGen_ZHbb_bb']==True))]
            #print(df)
            #
            df['ZH_pt'].clip(self.pt_bins[0]+1,self.pt_bins[-1]-1, inplace=True)
            df['pt_bin'] = pd.cut(df['ZH_pt'], bins=self.pt_bins, 
                                  labels=[f'Zhpt{i_bin}' for i_bin in range(len(self.pt_bins[:-1]))])
            # but we dont really care about pt_bin 0-200, so lets start at index 1
            self.hist_dict[y] = {}
            for i in range(1,len(self.pt_bins[:-1])):
                sub_df = df[df['pt_bin'] == f'Zhpt{i}']
                #quantiles = np.linspace(0,1,self.n_NN_bins+1) # actual quantiles, 10% intervals
                quantiles = [0.0, .05, .25, .35, .50, .70, 1.0] # actual quantiles, 10% intervals

                #if self.isblind:
                nn_df = sub_df[self.nn] 
                #print(nn_df)
                #else:
                #    nn_df = sub_df[self.nn][sub_df[self.nn]<=1.7]

                nn_bins = weighted_quantile(nn_df,
                                            quantiles, 
                                            getZhbbWeight(sub_df,y))
                nn_bins[0], nn_bins[-1] = 0,1 # explicity set bin edges to 0 and 1
                if self.isblind == False:
                    #nn_bins = nn_bins[:-1] # up to last bin but mass-side bands
                    nn_bins = nn_bins # dont do anything
                
                #nn_bins = nn_bins[1:] # drop first background dominated bin
                print(nn_bins)
                
                self.hist_dict[y][i] = functools.partial(
                    np.histogram2d,
                    bins=[nn_bins,(self.sdM_bins[i] if type(self.sdM_bins) is dict else self.sdM_bins)])
                #
            #
        #
                                            
    def __getitem__(self,y): # building this like a dictionary 
        try:
            return self.hist_dict[y]
        except KeyError:
            raise KeyError(f"{y} is not a valid top-level key!!!")