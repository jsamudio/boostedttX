import math
import awkward as ak
import numpy as np

'''
Helper to manually deal with weights and export in flat n-tuple
'''

def calc_weight(events, output, dataset, params, year):
    if events.metadata['isMC'] == 'True':
        print("LUMI:", params['lumi']['picobarns'][year]['tot'])
        norm_weight = (float(events.metadata['xsec'])*params['lumi']['picobarns'][year]['tot'])
        events['norm_weight'] = norm_weight
    else:
        events['norm_weight'] = 1

def add_weights_to_ttbb(events, sample): # Don't use this.
    # Hardcoded for now, what is the right way forward?
    if   '2L2Nu' in sample:
        events['norm_weight'] = 0.04799413040922868
    elif 'Semi' in sample:
        events['norm_weight'] = 0.10064983854787743
    elif 'Had' in sample:
        events['norm_weight'] = 0.14441797895972835

