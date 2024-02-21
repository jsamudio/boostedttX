import awkward as ak
import numpy as np

def zh_helper(events):
    deltaR = (lambda obj, obj2: ak.flatten(obj.metric_table(obj2), axis = 2))
    ZHCand = events.FatJetSorted[:,0]
    ak4 = events.JetGood
    '''
    Comparing ZH candidate with ak4, b, q ,l || l and b
    In general, want to output arrays of variables which are sorted by dR
    '''
    # ZH and AK4
    ZH_ak4_dr = ak.flatten(ZHCand.metric_table(ak4)) # This gives [ [ dR, ...], ...], innermost is comparing each obj1 with obj2
    # Argsort index array, placing nan first and then ascending by dR within 0.8
    ind_ZH_ak4_dr = ak.argsort(ak.where(ZH_ak4_dr <= 0.8, ZH_ak4_dr, np.nan), axis=1)
    # ak4 (in 0.8 dR to ZH cand.) btag sorted by dR, ascending with nan first
    ak4_btag_indRsort = ak.where(ZH_ak4_dr <= 0.8, ak4.btagDeepFlavB, np.nan)[ind_ZH_ak4_dr]




