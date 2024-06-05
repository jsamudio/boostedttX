import math
import awkward as ak
import numpy as np

'''
Helper to manually deal with weights and export in flat n-tuple
'''

def calc_weight(events, output, dataset, params):
    if events.metadata['isMC'] == 'True':
        norm_weight = (float(events.metadata['xsec'])*params.sample_params['lumi']['lumi']*1000)/output['sum_signOf_genweights'][dataset]
        events['norm_weight'] = norm_weight
    else:
        events['norm_weight'] = 1
