import math
import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

# We do this for various objects so why not make it a single function

deltaR = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_)))

def make_dRsorted_arr(obj1, obj2, comparison, proc):
    if (comparison == "gt"): #greaterthan
        nanArgsort = (lambda dr_: ak.argsort(ak.where(dr_ > 0.8, dr_, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ > 0.8, proc_, np.nan)[ind_])
    elif (comparison == "lt"):
        nanArgsort = (lambda dr_: ak.argsort(ak.where(dr_ <= 0.8, dr_, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ <= 0.8, proc_, np.nan)[ind_])
    else:
        nanArgsort = (lambda dr_: ak.argsort(dr_, axis=1))
        process_sorted = (lambda proc_, dr_, ind_: proc_[ind_])

    obj1_obj2_dr = deltaR(obj1, obj2)
    ind_obj1_obj2_dr = nanArgsort(obj1_obj2_dr)
    sorted_arr = ak.fill_none(ak.pad_none(process_sorted(proc, obj1_obj2_dr, ind_obj1_obj2_dr), 1, axis=1), np.nan)

    return sorted_arr

def make_ptsorted_arr(obj1, obj2, comparison, proc):
    if (comparison == "gt"): #greaterthan
        nanArgsort = (lambda dr_, obj1_: ak.argsort(-1*ak.where(dr_ > 0.8, obj1_.pt, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ > 0.8, proc_, np.nan)[ind_])
    else:
        nanArgsort = (lambda dr_, obj1_: ak.argsort(-1*ak.where(dr_ <= 0.8, obj1_.pt, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ <= 0.8, proc_, np.nan)[ind_])

    obj1_obj2_dr = deltaR(obj1, obj2)
    ind_obj1_obj2_pt = nanArgsort(obj1_obj2_dr, obj1)
    sorted_arr = ak.fill_none(ak.pad_none(process_sorted(proc, obj1_obj2_dr, ind_obj1_obj2_pt), 1, axis=1), np.nan)

    return sorted_arr

def zip_4vec(pt, eta, phi, m):
    vec = ak.zip(
            {
                "pt": pt,
                "eta": eta,
                "phi": phi,
                "mass": m,
            },
            with_name="PtEtaPhiMLorentzVector",
            behavior=vector.behavior,
    )
    return vec

def calc_SandA(pt_,eta_,phi_): # sphericity and aplanarity
    S_      = np.zeros((pt_.shape[0],3,3))
    pxyz_   = np.array([pt_*np.cos(phi_),pt_*np.sin(phi_),pt_*np.sinh(eta_)])
    p2_sum  = np.nansum(np.power(pt_*np.cosh(eta_),2),axis=1)
    for i_ in range(0,3):
        for j_ in range(0,3):
            S_[:,i_,j_] = np.nansum(pxyz_[i_]*pxyz_[j_], axis = 1)/p2_sum
    #
    eig_= -np.sort(-(np.linalg.eig(S_[:])[0]), axis = 1) # calc eig vals and sort by descending
    s_ = (3/2)*(eig_[:,1]+eig_[:,2])
    a_ = (3/2)*eig_[:,2]
    return s_, a_

def deltaPhi(phi1, phi2):
    try:
        dphi = np.subtract(phi1,phi2.T).T
    except (AttributeError) :
        dphi = phi1-phi2
    dphi_1 = ak.where(((dphi > math.pi) & (dphi != np.nan)), dphi - 2*math.pi, dphi)
    dphi_2 = ak.where(((dphi_1 <= -math.pi) & (dphi_1 != np.nan)), dphi_1 + 2*math.pi, dphi_1)

    return dphi_2

def calc_mtb(pt_, phi_, m_pt, m_phi):
    try:
        mtb2 =  2*(m_pt*pt_.T).T*(1 - np.cos(deltaPhi(m_phi,phi_)))
    except (AttributeError):
         mtb2 =  2*m_pt*pt_*(1 - np.cos(deltaPhi(phi_, m_phi)))
    return np.sqrt(mtb2)

def zh_helper(events):
    ZHCand = events.FatJetSorted[:,0]
    ak4 = events.JetGood
    ak8 = events.FatJetGood
    bjet = events.bJetGood
    qjet = events.qJetGood
    #lep = events.LeptonGood
    met = events.MET
    events['ZH_pt'] = ZHCand.pt
    events['ZH_M'] = ZHCand.particleNet_mass
    events['MET_pt'] = met.pt

    '''
    Dirty padding and array formatting to use old SandA calculation
    '''

    #spher, aplan = calc_SandA(
    #    np.append(ak.to_numpy(ak.fill_none(ak.pad_none(ak4.pt, max(ak.count(ak4.pt, axis = 1)), clip=True), np.nan)), ak.to_numpy(lep.pt), axis=1),
    #    np.append(ak.to_numpy(ak.fill_none(ak.pad_none(ak4.eta, max(ak.count(ak4.eta, axis = 1)), clip=True), np.nan)), ak.to_numpy(lep.eta), axis=1),
    #    np.append(ak.to_numpy(ak.fill_none(ak.pad_none(ak4.phi, max(ak.count(ak4.phi, axis = 1)), clip=True), np.nan)), ak.to_numpy(lep.phi), axis=1))

    #events["spher"] = spher
    #events["aplan"] = aplan

    '''
    Comparing ZH candidate with ak4, b, q ,l || l and b
    In general, want to output arrays of variables which are sorted by dR
    '''
    events["ZH_bbvLscore"] = ZHCand.particleNetMD_Xbb
    events["outZH_max_ak8pnetMass"] = ak.max(events.FatJetSorted.particleNet_mass[:,1:], axis=1)

    # ZH and AK4
    # ak4 (in 0.8 dR to ZH cand.) btag sorted by dR, ascending with nan first
    ak4_btag_indRsort = make_dRsorted_arr(ZHCand, ak4, "lt", ak4.btagDeepFlavB)

    events["ak4_bestb_inZH"] = ak.max(ak4_btag_indRsort, axis=1)
    events["ak4_worstb_inZH"] = ak.min(ak4_btag_indRsort, axis=1)

    # ZH and b

    b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort = [make_dRsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_dRsort_vec = zip_4vec(b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort)

    ind_close_b = ak.argmax(ak.nan_to_num(b_pt_dRsort, nan=-1), axis = 1, keepdims=True)

    close_b = b_dRsort_vec[ind_close_b]

    events["ZH_closeb_invM"] = ak.flatten((ZHCand + close_b).mass)

    b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort, b_btag_ptsort = [make_ptsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass, bjet.btagDeepFlavB]]
    b_ptsort_vec = zip_4vec(b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort)
    #print(ak.pad_none(b_ptsort_vec, 1, axis=1))

    # ZH and q

    q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort = [make_dRsorted_arr(ZHCand, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass]]
    q_dRsort_vec = zip_4vec(q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort)

    q_pt_ptsort, q_eta_ptsort, q_phi_ptsort, q_mass_ptsort, q_btag_ptsort = [make_ptsorted_arr(ZHCand, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass, qjet.btagDeepFlavB]]
    q_ptsort_vec = zip_4vec(q_pt_ptsort, q_eta_ptsort, q_phi_ptsort, q_mass_ptsort)

    # Combinations

    #ZH_l_dr = deltaR(ZHCand, lep)

    #ZH_l_invM = (ZHCand + lep).mass #might need this to be specifically pNet mass

    #events["ZH_l_dr"] = ak.flatten(ZH_l_dr)
    #events["ZH_l_invM"] = ak.flatten(ZH_l_invM)

    # Leading-Subleading (pt) b and q combinations

    ind_leading_b = ak.argmax(ak.nan_to_num(b_pt_ptsort, nan=-1), axis = 1, keepdims=True)

    leading_b = b_ptsort_vec[ind_leading_b]

    leading_b_btag = b_btag_ptsort[ind_leading_b]

    ind_subLeading_b = np.where(ind_leading_b + 1 < ak.count(b_pt_ptsort, axis=1), ind_leading_b + 1, ind_leading_b)

    subLeading_b = b_ptsort_vec[ind_subLeading_b]

    subLeading_b_btag = b_btag_ptsort[ind_subLeading_b]

    b_b_dr = deltaR(leading_b, subLeading_b)

    ind_leading_q = ak.argmax(ak.nan_to_num(q_pt_ptsort, nan=-1), axis = 1, keepdims=True)

    leading_q = q_ptsort_vec[ind_leading_q]

    leading_q_btag = q_btag_ptsort[ind_leading_q]

    ind_subLeading_q = np.where(ind_leading_q + 1 < ak.count(q_pt_ptsort, axis=1), ind_leading_q + 1, ind_leading_q)

    subLeading_q = q_ptsort_vec[ind_subLeading_q]

    subLeading_q_btag = q_btag_ptsort[ind_subLeading_q]

    q_q_dr = deltaR(leading_q, subLeading_q)

    # Nearest and second nearest b to l
    # Any NaN value is a b within the dR cone of the ZH candidate

    b_pt_dRsort_l, b_eta_dRsort_l, b_phi_dRsort_l, b_mass_dRsort_l = [make_dRsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_dRsort_l_vec = zip_4vec(b_pt_dRsort_l, b_eta_dRsort_l, b_phi_dRsort_l, b_mass_dRsort_l)

    ind_near_b = ak.argmax(ak.nan_to_num(b_pt_dRsort_l, nan=-1), axis = 1, keepdims=True)

    near_b = b_dRsort_l_vec[ind_near_b]

    ind_subNear_b = np.where(ind_near_b + 1 < ak.count(b_pt_dRsort_l, axis=1), ind_near_b + 1, ind_near_b)

    subNear_b = b_dRsort_l_vec[ind_subNear_b]

    #l_b1_invM = (lep + near_b).mass
    #l_b1_dr = deltaR(lep, near_b)

    #l_b2_invM = (lep + subNear_b).mass
    #l_b2_dr = deltaR(lep, subNear_b)

    #l_b1_mtb = calc_mtb((lep+leading_b).pt, (lep+leading_b).phi, met.pt, met.phi)
    #l_b2_mtb = calc_mtb((lep+leading_b).pt, (lep+leading_b).phi, met.pt, met.phi)

    #events["l_b1_invM"] = ak.flatten(l_b1_invM)
    #events["l_b2_invM"] = ak.flatten(l_b2_invM)

    #events["l_b1_dr"] = ak.flatten(l_b1_dr)
    #events["l_b2_dr"] = ak.flatten(l_b2_dr)

    #events["l_b2_mtb"] = l_b2_mtb
    # Combinations involving ZH and leading(subleading) b and q

    events["outZH_b1_pt"] = ak.flatten(leading_b.pt)
    events["outZH_b2_pt"] = ak.flatten(subLeading_b.pt)
    events["outZH_b1_score"] = ak.flatten(leading_b_btag)
    events["outZH_b2_score"] = ak.flatten(subLeading_b_btag)
    events["outZH_q1_pt"] = ak.flatten(leading_q.pt)
    events["outZH_q2_pt"] = ak.flatten(subLeading_q.pt)
    events["outZH_q1_score"] = ak.flatten(leading_q_btag)
    events["outZH_q2_score"] = ak.flatten(subLeading_q_btag)

    b1_q_dr = deltaR(leading_b, q_ptsort_vec)
    b2_q_dr = deltaR(subLeading_b, q_ptsort_vec)
    # for now, nanmin not implemented so events with only nan, place inf
    events["outZH_b1_q_mindr"] = ak.min(ak.where(b1_q_dr == np.nan, np.inf, b1_q_dr), axis=1)
    events["outZH_b2_q_mindr"] = ak.min(ak.where(b2_q_dr == np.nan, np.inf, b2_q_dr), axis=1)

    #events["l_b1_mtb"] = ak.flatten(l_b1_mtb)
    #events["l_b2_mtb"] = ak.flatten(l_b2_mtb)

    # q near b
    #print(leading_b)
    #print(qjet)
    q_pt_b1dRsort, q_eta_b1dRsort, q_phi_b1dRsort, q_mass_b1dRsort = [make_dRsorted_arr(leading_b, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass]]
    q_b1dRsort_vec = zip_4vec(q_pt_b1dRsort, q_eta_b1dRsort, q_phi_b1dRsort, q_mass_b1dRsort)

    q_pt_b2dRsort, q_eta_b2dRsort, q_phi_b2dRsort, q_mass_b2dRsort = [make_dRsorted_arr(subLeading_b, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass]]
    q_b2dRsort_vec = zip_4vec(q_pt_b2dRsort, q_eta_b2dRsort, q_phi_b2dRsort, q_mass_b2dRsort)
    # nearest and second nearest

    ind_nearb1_q = ak.argmax(ak.nan_to_num(q_pt_b1dRsort, nan=-1), axis = 1, keepdims=True)

    nearb1_q = q_b1dRsort_vec[ind_nearb1_q]

    ind_subNearb1_q = np.where(ind_nearb1_q + 1 < ak.count(q_pt_b1dRsort, axis=1), ind_nearb1_q + 1, ind_nearb1_q)


    subNearb1_q = q_b1dRsort_vec[ind_subNearb1_q]

    events["outZH_q_q_dr_nearb1"] = ak.flatten(deltaR(nearb1_q, subNearb1_q))

    ind_nearb2_q = ak.argmax(ak.nan_to_num(q_pt_b2dRsort, nan=-1), axis = 1, keepdims=True)

    nearb2_q = q_b2dRsort_vec[ind_nearb2_q]

    ind_subNearb2_q = np.where(ind_nearb2_q + 1 < ak.count(q_pt_b2dRsort, axis=1), ind_nearb2_q + 1, ind_nearb2_q)

    subNearb2_q = q_b2dRsort_vec[ind_subNearb2_q]

    events["outZH_q_q_dr_nearb2"] = ak.flatten(deltaR(nearb2_q, subNearb2_q))

    events["outZH_qq_M_nearb1"] = ak.flatten((nearb1_q + subNearb1_q).mass)
    events["outZH_qq_M_nearb2"] = ak.flatten((nearb2_q + subNearb2_q).mass)
    events["outZH_b1q_M"] = ak.flatten((leading_b + nearb1_q).mass)
    events["outZH_b1q_M"] = ak.flatten((subLeading_b + nearb2_q).mass)

    b1qq = (leading_b + nearb1_q + subNearb1_q)
    b2qq = (subLeading_b + nearb2_q + subNearb2_q)

    events["outZH_b1_qq_dr"] = ak.flatten(deltaR(leading_b, (nearb1_q + subNearb1_q)))
    events["outZH_b2_qq_dr"] = ak.flatten(deltaR(subLeading_b, (nearb2_q + subNearb2_q)))
    events["outZH_b1qq_M"] = ak.flatten(b1qq.mass)
    events["outZH_b2qq_M"] = ak.flatten(b2qq.mass)
    events["ZH_b1qq_dr"] = ak.flatten(deltaR(ZHCand, b1qq))
    events["ZH_b2qq_dr"] = ak.flatten(deltaR(ZHCand, b2qq))

    #lbb1qq = (lep + leading_b + subLeading_b + nearb1_q + subNearb1_q)
    #lbb2qq = (lep + leading_b + subLeading_b + nearb2_q + subNearb2_q)

    #events["ZH_lbb1qq_dr"] = ak.flatten(deltaR(ZHCand, lbb1qq))
    #events["ZH_lbb2qq_dr"] = ak.flatten(deltaR(ZHCand, lbb2qq))

    ZH_b_dr = deltaR(ZHCand, bjet)
    ZH_q_dr = deltaR(ZHCand, qjet)
    n_b_outZhbb = ak.sum(ZH_b_dr > 0.8, axis = 1)
    n_b_inZhbb = ak.sum(ZH_b_dr <= 0.8, axis = 1)

    n_q_outZhbb = ak.sum(ZH_q_dr > 0.8, axis = 1)
    n_q_inZhbb = ak.sum(ZH_q_dr <= 0.8, axis = 1)

    events["n_b_inZH"] = n_b_inZhbb
    events["n_b_outZH"] = n_b_outZhbb
    events["n_q_inZH"] = n_q_inZhbb
    events["n_q_outZH"] = n_q_outZhbb

    b_pt_indRsort, b_eta_indRsort, b_phi_indRsort, b_mass_indRsort = [make_dRsorted_arr(ZHCand, bjet, "lt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_indRsort_vec = zip_4vec(b_pt_indRsort, b_eta_indRsort, b_phi_indRsort, b_mass_indRsort)

    ind_b_outZH = ak.argmax(ak.nan_to_num(b_pt_dRsort, nan=-1), axis = 1, keepdims=True)

    outZH_b = b_dRsort_vec[ind_b_outZH]

    ind_b_inZH = ak.argmax(ak.nan_to_num(b_pt_indRsort, nan=-1), axis = 1, keepdims=True)

    inZH_b = b_indRsort_vec[ind_b_inZH]

    events["inZHb_outZHb_dr"] = ak.flatten(deltaR(inZH_b, outZH_b))

    ht_b = ak.sum(ak.nan_to_num(b_pt_ptsort, nan=0), axis=1)
    #sc_pt_outZH = ht_b + ak.sum(ak.nan_to_num(q_pt_ptsort, nan=0), axis=1) + lep.pt
    events["ht_b"] = ht_b
    #events["ht_outZH"] = ak.flatten(sc_pt_outZH)

    events["outZH_b12_m"] = ak.flatten((leading_b + subLeading_b).mass)
    events["outZH_b12_dr"] = ak.flatten(deltaR(leading_b, subLeading_b))
    events["nonZHbb_q1_dr"] = ak.flatten(deltaR(ZHCand, leading_q))
    events["nonZHbb_b1_dr"] = ak.flatten(deltaR(ZHCand, leading_b))

    events["n_ak4jets"] = ak.count(ak4.pt, axis=1)
    events["n_ak8jets"] = ak.count(ak8.pt, axis=1)
    events["n_ak8_ZHbb"] = ak.sum(ak.where(ak8.particleNetMD_Xbb > 0.8, 1, 0), axis=1)
