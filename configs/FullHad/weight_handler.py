import math
import awkward as ak
import numpy as np

'''
Helper to manually deal with weights and export in flat n-tuple
'''

def calc_weight(events, output, dataset, params):
    if events.metadata['isMC'] == 'True':
        norm_weight = (float(events.metadata['xsec'])*params.sample_params['lumi']['lumi']*1000)
        events['norm_weight'] = norm_weight
    else:
        events['norm_weight'] = 1
# Need to manually set weights for TTbb samples. These come from process norms in the old framework, but will need application here.
# Also need to apply the toppt reweighting as well because this may have the effect that we are looking for

def add_weights_to_ttbb(events, sample):
    # Hardcoded for now, what is the right way forward?
    if   '2L2Nu' in sample:
        events['norm_weight'] = 0.04799413040922868
    elif 'Semi' in sample:
        events['norm_weight'] = 0.10064983854787743
    elif 'Had' in sample:
        events['norm_weight'] = 0.14441797895972835

