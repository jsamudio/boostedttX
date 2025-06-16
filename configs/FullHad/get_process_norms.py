import sys
import numpy as np
import pandas as pd
import json
import re
from glob import glob

'''
This needs to accomplish a few things, first and foremost calcuating the total yield from each category and writing to a json file.
Systematics are needed too but we can treat those once we have them. For now we should focus on saving the following:
- weight
- yield
- stxs yield
- xsec

Format should be:
{
    "Year": {
        "Process": {
            "yield": xyz,
        },
    },
}

With signal as pt bins ttH0/1/2/3 and inclusive ttH.
Backgrounds are lumped together by process except for ttbb;
add an additional category per subsample with weight and xsec
'''

processes = ['TTBar', 'tt_B', 'ttZ', 'ttH']
content = ['yield']
pt_bins = [0,200,300,450]

def main():
    json_dict = {}
    df = read_pkl()
    # for loop covers the year entry
    for y in ['2017']:
        # worker should be replaced with the proper data entry
        json_dict[y] = __worker(df, y)
    out_json_file = './process_norms/process_norms_ttbbw_run2.json'
    with open(out_json_file, 'w') as jsf:
        json.dump(json_dict, jsf, indent=4)

def read_pkl(): 
    get_pickle = (lambda s: pd.read_pickle(f'pickled/SpanetInferenceDouble_{s}.pkl'))
    df = pd.concat([get_pickle(s) for s in processes], axis='rows', ignore_index=True)
    df['pt_bin'] = pd.cut(df['ZH_pt'], bins=pt_bins+[np.inf], labels=[i_bin for i_bin in range(len(pt_bins))])
    return df

# return array of weights; will need subprocess split for ttbb
def get_weight(process):
    return df[df['process'] == process]['norm_weight'].as_array(dtype='float')

# should be norm_weight * sign of genweight for the total events in the sample, not just post-cut
def get_tot_weight(df_):
    tot_weight = (df_['norm_weight'] * np.sign(df_['genWeight']) * df_['topptWeight']) # and other weights
    return tot_weight

def get_yield(df_, process):
    df = df_[df_['process'] == process]
    tot_weight = get_tot_weight(df)
    event_yield = np.sum(tot_weight)
    return event_yield

def divide_signal(df_, process, content_dict, process_dict):
    if process in ['ttZ', 'ttH']:
        content_dict = {}
        for i,_ in enumerate(pt_bins[:-1]):
            new_sig_name = f"{process}{i}"
            content_dict['yield'] = get_yield(df_[(df_['ZH_pt']>= pt_bins[i]) &
                                          (df_['ZH_pt'] < pt_bins[i+1])], process)
            process_dict[new_sig_name] = content_dict
            content_dict = {}
        content_dict['yield'] = get_yield(df_[df_['ZH_pt'] > pt_bins[-1]], process)
        process_dict[f'{process}{len(pt_bins)-1}'] = content_dict
    return process_dict

def __worker(df_, y):
    process_dict = {}
    for i in processes:
        content_dict = {}
        content_dict['yield'] = get_yield(df_, i)
        #content_dict['otherVar'] = "test"
        process_dict[f'{i}'] = content_dict
        divide_signal(df_, i, content_dict, process_dict)
        print(process_dict)
    return process_dict

if __name__=='__main__':
    main()