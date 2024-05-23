import math
import awkward as ak
import numpy as np

'''
Helper to manually deal with weights and export in flat n-tuple
'''

def calc_weight(events, sample_params):
    if events.metadata['isMC'] == 'True':
        norm_weight = (events.metadata['xsec']*sample_params['kf']*sample_params['lumi']*1000)/events['sum_sign_genw']
        events['norm_weight'] = norm_weight
    else
        events['norm_weight'] = 1
