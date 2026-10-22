import math
import awkward as ak
import numpy as np
import pandas as pd
from coffea.nanoevents.methods import vector, nanoaod

deltaR = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_)))

def match_gen_lep(events):
    gen_id = events.GenPart.pdgId
    gen_mom = events.GenPart.genPartIdxMother
    gen_st = events.GenPart.status
    gen_pt = events.GenPart.pt
    #islep = (((abs(gen_id) == 11) | (abs(gen_id) == 13)) & ((abs(gen_id[gen_mom[gen_mom]]) == 6) | (abs(gen_id[gen_mom[gen_mom]]) == 24)) &(abs(gen_id[gen_mom]) ==24))
    islep = (((abs(gen_id) == 11) | (abs(gen_id) == 13)) & (abs(events.GenPart[gen_mom].distinctParent.pdgId) == 6) &(abs(gen_id[gen_mom]) ==24))
    lep_match_dr = deltaR(events.LeptonGood, events.GenPart[islep])
    events["matchedGenLep"] = (ak.sum(lep_match_dr <= 0.1, axis=1) > 0)

def match_gen_tt(events, sample):
    # first get tt type
    if   '2L2Nu' in sample:
        events['tt_type'] = 'Di'
    elif 'Semi' in sample:
        events['tt_type'] = 'Semi'
    elif 'Had' in sample:
        events['tt_type'] = 'Had'

    gen_id = events.GenPart.pdgId
    gen_mom = events.GenPart.genPartIdxMother
    gen_pt = events.GenPart.pt

    gentt_bb = events.genTtbarId
    gentt_bb = gentt_bb % 100
    is_tt_C = ( (gentt_bb>=41) & (gentt_bb<50) )
    is_tt_B = gentt_bb>=51

    is_tt_b  = gentt_bb == 51
    is_tt_2b = gentt_bb == 52
    is_tt_bb = gentt_bb >= 53
    #
    is_tt_c  = gentt_bb == 41
    is_tt_2c = gentt_bb == 42
    is_tt_cc = (gentt_bb >= 43) & (gentt_bb<50)
    #
    events['tt_C'] = is_tt_C
    events['tt_B'] = is_tt_B
    events['tt_b'] = is_tt_b
    events['tt_2b'] = is_tt_2b
    events['tt_bb'] = is_tt_bb
    #
    events['tt_c'] = is_tt_c
    events['tt_2c'] = is_tt_2c
    events['tt_cc'] = is_tt_cc
    if 'TTTo' in sample:
        events['process'] = np.where(events['tt_B'] == True, 'old_tt_B', 'TTBar')
    elif (('TTbb' in sample) | ('ttbb' in sample)):
        events['process'] = np.where(events['tt_B'] == True, 'tt_B', 'non_tt_B')

    # check if a Z is from the hard process
    events['has_Z'] = ak.sum((gen_mom <= 0) & (abs(gen_id) == 23), axis=1) > 0
    events['has_H'] = ak.sum((gen_mom <= 0) & (abs(gen_id) == 25), axis=1) > 0

    # calculate toppt weight for powheg only
    # last-copy hard-process top quarks (status 62)
    tt_pt = gen_pt[(abs(gen_id) == 6) & (events.GenPart.status == 62)]
    # Using the newer theo (NNLO QCD + NLO EW) corrections which is better for BSM analysis aspects
    sf = (lambda x: 0.103*np.exp(-0.0118*np.clip(x,0,np.inf)) - 0.000134*np.clip(x,0,np.inf) + 0.973) #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_3_The_Effective_Field_The
    # extrapolation from 13 to 13.6 TeV (TOP-25-018), multiplying the 13 TeV correction
    sf13p6 = (lambda x: 0.985 + 0.0004*np.clip(x,0,np.inf))
    #https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf'
    # the theo toppt event re-weighting unc. is based on [1, w**2] where w is the event reweighting
    # per-top weight SF13(pT) * SF13.6(pT); event weight = geometric mean over top and antitop
    per_top = (lambda x: sf(x) * sf13p6(x))
    toppt_rwgt = np.sqrt(per_top(tt_pt[:,0]) * per_top(tt_pt[:,1]))
    toppt_rwgt_up = np.where(toppt_rwgt > 1.0, toppt_rwgt**2,  1.0)
    toppt_rwgt_dn = np.where(toppt_rwgt < 1.0, toppt_rwgt**2,  1.0)
    events['tt_pt1'], events['tt_pt2'] = tt_pt[:,0], tt_pt[:,1]
    events['topptWeight']      = toppt_rwgt
    events['topptWeight_Up']   = toppt_rwgt_up
    events['topptWeight_Down'] = toppt_rwgt_dn

def fill1e(one_d_array):
    return pd.DataFrame.from_records(one_d_array).values.flatten()

def fillne(nd_array):
    counts_ = np.array(ak.count(nd_array))
    rect_ , flat_ = np.full((len(nd_array),max(counts_)),np.nan), np.array(ak.flatten(nd_array))
    def cfillne(ja, o, c):
        rows, _ = o.shape
        c_i = 0
        for i in range(rows):
            o[i,:c[i]] = ja[c_i:c_i+c[i]]
            c_i += c[i]
        return o
    return cfillne(flat_, rect_, counts_)
    #

def match_gen_sig(events, sample):
    ZHCand = events.FatJetSorted[:,0]
    gen_id = events.GenPart.pdgId
    gen_mom = events.GenPart.genPartIdxMother
    gen_pt = events.GenPart.pt
    gen_st = events.GenPart.status

    istt = (abs(gen_id) == 6)
    #
    isbb_fromZ     = ((abs(gen_id) == 5) & (gen_id[gen_mom] == 23) & (gen_st[gen_mom] == 62))
    isqq_fromZ     = ((abs(gen_id) <  5) & (gen_id[gen_mom] == 23) & (gen_st[gen_mom] == 62))
    isllnunu_fromZ = ((abs(gen_id) >=  11) & (abs(gen_id) <= 16) & (gen_id[gen_mom] == 23) & (gen_st[gen_mom] == 62))
    isbb_fromH    = ((abs(gen_id) == 5) & (gen_id[gen_mom] == 25) & (gen_st[gen_mom] == 62))
    isnonbb_fromH = ((abs(gen_id) != 5) & (gen_id[gen_mom] == 25) & (gen_st[gen_mom] == 62))
    isHbb     = ((gen_id == 25) & (ak.sum(isbb_fromH, axis=1) == 2) & (gen_st == 62))
    isHnonbb  = ((gen_id == 25) & (ak.sum(isbb_fromH, axis=1) == 0) & (gen_st == 62))
    isZbb  = ((gen_id == 23) & (ak.sum(isbb_fromZ, axis=1) == 2) & ((ak.sum(isHbb, axis=1) == 0) & (ak.sum(isHnonbb, axis=1) == 0)) & (gen_st == 62))
    isZqq  = ((gen_id == 23) & (ak.sum(isqq_fromZ, axis=1) == 2) & ((ak.sum(isHbb, axis=1) == 0) & (ak.sum(isHnonbb, axis=1) == 0)) & (gen_st == 62))
    isZllnunu = ((gen_id == 23) & (ak.sum(isllnunu_fromZ, axis=1) == 2) & ((ak.sum(isHbb, axis=1) == 0) & (ak.sum(isHnonbb, axis=1) == 0)) & (gen_st == 62))
    isZH = ((isHbb) | (isZbb) | (isZqq) | (isZllnunu) | (isHnonbb))
    #print("recoZH",len(rZh_eta),len(rZh_phi))

    zh_match_dR = deltaR(ZHCand, events.GenPart[(isZH)])
    rzh_matchb_dR = deltaR(ZHCand, events.GenPart[(isbb_fromZ) | (isbb_fromH)])
    rzh_matchtt_dR = deltaR(ZHCand, events.GenPart[(istt)])
    zh_matchbb = (ak.sum(rzh_matchb_dR <= 0.6, axis=1) == 2)
    zh_matchb = (ak.sum(rzh_matchb_dR <= 0.6, axis=1) == 1)
    zh_nomatchb = (ak.sum(rzh_matchb_dR <= 0.6, axis=1) == 0)

    zh_match = ((zh_match_dR <= 0.6) & (events.GenPart[isZH].pt >= 100) & (abs(events.GenPart[isZH].eta) <= 2.4))

    events['Zbb']= (ak.sum(isZbb, axis=1) > 0)
    events['Hbb']= (ak.sum(isHbb, axis=1) > 0)
    events['Hnonbb']= (ak.sum(isHnonbb, axis=1) > 0)
    events['Zqq']= (ak.sum(isZqq, axis=1) > 0)
    events['Zllnunu']= (ak.sum(isZllnunu, axis=1) > 0)

    #zh_pt = fill1e(gen_pt[isZH]).flatten()
    events['genZHpt'] = ak.fill_none(ak.firsts(gen_pt[isZH]), -1.0)
    #events['genZHpt'] = zh_pt

    events['matchedGenZH']    = ak.sum(zh_match, axis=1) > 0
    events['matchedGen_Zbb']  = ((ak.sum(zh_match, axis=1) > 0) & (events['matchedGenLep']) & (ak.sum(isZbb,axis=1) >  0))
    events['matchedGen_Hbb']  = ((ak.sum(zh_match, axis=1) > 0) & (events['matchedGenLep']) & (ak.sum(isHbb,axis=1) >  0))
    #print(sum(events['matchedGen_Hbb']))
    #events['matchedGen_Hbb2']  = ((ak.sum(zh_match, axis=1) > 0) & (events['matchedGenLep2']) & (ak.sum(isHbb,axis=1) >  0))
    #print(sum(events['matchedGen_Hbb2']))
    events['matchedGen_ZHbb'] = ((ak.sum(zh_match, axis=1) > 0) & (events['matchedGenLep']) & ((ak.sum(isZbb,axis=1) + ak.sum(isHbb, axis=1)) > 0))
    events['matchedGen_Zqq']  = ((ak.sum(zh_match, axis=1) > 0) & (events['matchedGenLep']) & (ak.sum(isZqq,axis=1) >  0))
    #
    events['matchedGen_ZHbb_bb']  = ((events['matchedGen_ZHbb'] == 1) & (zh_matchbb  == 1))
    events['matchedGen_ZHbb_b']   = ((events['matchedGen_ZHbb'] == 1) & (zh_matchb   == 1))
    events['matchedGen_ZHbb_nob'] = ((events['matchedGen_ZHbb'] == 1) & (zh_nomatchb == 1))

    # Run 3: the process label comes from samples.PROCESS_LABELS (set by the workflow). There is
    # no dedicated TTZ->bb sample, so TTZToQQ Z->bb events stay in the ttZ signal.

    events['topptWeight']      = 1
    events['topptWeight_Down'] = 1
    events['topptWeight_Up'] = 1

def match_gen_bb_to_zh(events):
    """
    Evaluates geometric gen-matching for DeepAK8 MD Xbb scale factors.
    Flags events as 'real' (True) if exactly 2 gen b-quarks are within dR < 0.6 
    of the Z/H candidate, otherwise 'fake' (False).
    """
    # Assuming the first sorted fat jet is your Z/H candidate
    ZHCand = events.FatJetSorted[:, 0]
    
    # Isolate gen b-quarks
    gen_id = events.GenPart.pdgId
    
    # Safely get the parent pdgId using coffea NanoEvents behavior. 
    # Fill none with 0 for initial partons that have no parent.
    parent_id = ak.fill_none(events.GenPart.parent.pdgId, 0)
    
    # Requirement: Particle is a b-quark, and its parent is NOT a b-quark 
    # (prevents double counting intermediate shower copies)
    is_b = (abs(gen_id) == 5) & (abs(parent_id) != 5)
    
    # Calculate dR using the provided lambda function
    b_match_dr = deltaR(ZHCand, events.GenPart[is_b])
    
    # Count how many gen b-quarks fall within the 0.6 dR cone
    events['n_genb_matchZH'] = ak.sum(b_match_dr <= 0.6, axis=-1)
    
    # Apply the strict BTV definition for the scale factor phase space:
    # 'Real' = Exactly 2 b-quarks. 'Fake' = Anything else.
    events['xbb_genmatch'] = (events['n_genb_matchZH'] == 2)

def match_tt_products(events):
    ZHCand = events.FatJetSorted[:,0]
    gen_id = events.GenPart.pdgId
    gen_mom = events.GenPart.genPartIdxMother
    gen_pt = events.GenPart.pt
    gen_st = events.GenPart.status
    
    # Look at other decay products dR with respect to the ZH candidate

    ispart_fromWfromtop = (((abs(gen_id) < 5)) & 
                                   ((abs(gen_id[gen_mom[gen_mom]]) == 6) & 
                                    (abs(gen_id[gen_mom]) ==24)))
    
    isbb_fromtt     = ((abs(gen_id) == 5) & (abs(gen_id[gen_mom]) == 6))

    events['min_wpart_ZH_dR'] = ak.min(deltaR(ZHCand, events.GenPart[(ispart_fromWfromtop)]), axis=-1)
    events['max_wpart_ZH_dR'] = ak.max(deltaR(ZHCand, events.GenPart[(ispart_fromWfromtop)]), axis=-1)
    events['n_wpart_ZH_dR_0p4'] = ak.sum(deltaR(ZHCand, events.GenPart[(ispart_fromWfromtop)]) > 0.4, axis=-1)
    events['n_wpart_ZH_dR_0p8'] = ak.sum(deltaR(ZHCand, events.GenPart[(ispart_fromWfromtop)]) > 0.8, axis=-1)
    events['n_wpart_ZH_dR_1p2'] = ak.sum(deltaR(ZHCand, events.GenPart[(ispart_fromWfromtop)]) > 1.2, axis=-1)
    #print("Delta R: ", events.wpart_ZH_dR)
    events['min_topb_ZH_dR'] = ak.min(deltaR(ZHCand, events.GenPart[(isbb_fromtt)]), axis=-1)
    events['max_topb_ZH_dR'] = ak.max(deltaR(ZHCand, events.GenPart[(isbb_fromtt)]), axis=-1)
    events['n_topb_ZH_dR_0p4'] = ak.sum(deltaR(ZHCand, events.GenPart[(isbb_fromtt)]) > 0.4, axis=-1)
    events['n_topb_ZH_dR_0p8'] = ak.sum(deltaR(ZHCand, events.GenPart[(isbb_fromtt)]) > 0.8, axis=-1)
    events['n_topb_ZH_dR_1p2'] = ak.sum(deltaR(ZHCand, events.GenPart[(isbb_fromtt)]) > 1.2, axis=-1)
    
