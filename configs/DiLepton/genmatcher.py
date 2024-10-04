import math
import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector, nanoaod

deltaR = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_, axis=1)))

def match_gen_lep(events):
    gen_id = events.GenPart.pdgId
    gen_mom = events.GenPart.genPartIdxMother
    gen_st = events.GenPart.status
    gen_pt = events.GenPart.pt
    #islep = (((abs(gen_id) == 11) | (abs(gen_id) == 13)) & ((abs(gen_id[gen_mom[gen_mom]]) == 6) | (abs(gen_id[gen_mom[gen_mom]]) == 24)) &(abs(gen_id[gen_mom]) ==24))
    islep = (((abs(gen_id) == 11) | (abs(gen_id) == 13)) & (abs(events.GenPart[gen_mom].distinctParent.pdgId) == 6) &(abs(gen_id[gen_mom]) ==24))
    #deltaR2 = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_, axis=None)))
    #lep_match_dr = deltaR2(events.LeptonGoodDi, events.GenPart[islep])
    print(len(events.LeptonGoodDi))
    print(len(events.GenPart[islep]))
    lep_match_dr = ak.flatten(events.LeptonGoodDi.metric_table(events.GenPart[islep]), axis = -1)
    print(lep_match_dr)
    print(len(ak.sum(lep_match_dr <= 0.1, axis=-1) > 1))
    print(sum(ak.sum(lep_match_dr <= 0.1, axis=-1) > 1))
    print(len(events.GenPart))
    events["matchedGenLep"] = (ak.sum(lep_match_dr <= 0.1, axis=-1) > 1)
    #events["matchedGenLep"] = True

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
    elif 'TTbb' in sample:
        events['process'] = np.where(events['tt_B'] == True, 'tt_B', 'non_tt_B')
    print(events['process'])
    print('Total',len(is_tt_B))
    print('tt+B', sum(is_tt_B))
    print('tt+b', sum(is_tt_b))
    print('tt+2b', sum(is_tt_2b))
    print('tt+bb', sum(is_tt_bb))

    # check if a Z is from the hard process
    events['has_Z'] = ak.sum((gen_mom <= 0) & (abs(gen_id) == 23), axis=1) > 0
    events['has_H'] = ak.sum((gen_mom <= 0) & (abs(gen_id) == 25), axis=1) > 0
    print('has Z, <0 mom ',sum(events['has_Z']))
    print('has H, <0 mom ',sum(events['has_H']))
    print('has Z ',sum(ak.sum(abs(gen_id) == 23, axis=1) >0))
    print('has H ',sum(ak.sum(abs(gen_id) == 25, axis=1) >0))
    print('No Z Total', len(is_tt_B [events['has_Z'] == False]))
    print('No Z tt+B',  sum(is_tt_B [events['has_Z'] == False]))
    print('No Z tt+b',  sum(is_tt_b [events['has_Z'] == False]))
    print('No Z tt+2b', sum(is_tt_2b[events['has_Z'] == False]))
    print('No Z tt+bb', sum(is_tt_bb[events['has_Z'] == False]))

    # calculate toppt weight for powheg only
    tt_pt = gen_pt[(abs(gen_id) == 6)]
    # Using the newer theo (NNLO QCD + NLO EW) corrections which is better for BSM analysis aspects
    sf = (lambda x: 0.103*np.exp(-0.0118*np.clip(x,0,np.inf)) - 0.000134*np.clip(x,0,np.inf) + 0.973) #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Case_3_3_The_Effective_Field_The

    #https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf'
    # the theo toppt event re-weighting unc. is based on [1, w**2] where w is the event reweighting
    toppt_rwgt = np.sqrt(sf(tt_pt[:,0]) * sf(tt_pt[:,1]))
    toppt_rwgt_up = np.where(toppt_rwgt > 1.0, toppt_rwgt**2,  1.0)
    toppt_rwgt_dn = np.where(toppt_rwgt < 1.0, toppt_rwgt**2,  1.0)
    events['tt_pt1'], events['tt_pt2'] = tt_pt[:,0], tt_pt[:,1]
    events['topptWeight']      = toppt_rwgt
    events['topptWeight_Up']   = toppt_rwgt_up
    events['topptWeight_Down'] = toppt_rwgt_dn

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
    print("isHbb", sum(ak.sum(isHbb, axis=1)))
    print("isHnonbb", sum(ak.sum(isHnonbb, axis=1)))
    print("isZbb axis 1", sum(ak.sum(isZbb,axis=1)))
    print("isZbb axis -1", sum(ak.sum(isZbb,axis=-1)))
    print("isZqq", sum(ak.sum(isZqq, axis=1)))
    print("isZllnunu", sum(ak.sum(isZllnunu, axis=1)))
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

    if 'TTZToQQ' in sample:
        events['process'] = np.where(events['Zbb'] == True, 'old_ttZbb', 'ttZ')
    elif 'TTZ' in sample:
        events['process'] = 'ttZ'
    elif 'ttH' in sample:
        events['process'] = 'ttH'

    events['topptWeight']      = 1
