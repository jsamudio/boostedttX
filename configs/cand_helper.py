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
    sorted_arr = process_sorted(proc, obj1_obj2_dr, ind_obj1_obj2_dr)

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
    sorted_arr = process_sorted(proc, obj1_obj2_dr, ind_obj1_obj2_pt)

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

def invM(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2):
    m1sq = np.power(m1,2)
    m2sq = np.power(m2,2)
    pt1cosheta1_sq = np.power(pt1,2)*np.power(np.cosh(eta1),2)
    pt2cosheta2_sq = np.power(pt2,2)*np.power(np.cosh(eta2),2)
    try:
        E1pE22 = np.power((np.sqrt(pt1cosheta1_sq+m1sq)+np.sqrt(pt2cosheta2_sq+m2sq).T).T,2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = (np.sinh(eta1)*np.sinh(eta2).T).T
        p1dotp2 = (pt1*pt2.T).T*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - (pt1cosheta1_sq + pt2cosheta2_sq.T).T - 2*p1dotp2
    except (AttributeError):
        E1pE22 = np.power(np.sqrt(pt1cosheta1_sq+m1sq)+np.sqrt(pt2cosheta2_sq+m2sq),2)
        cosphi1phi2 = np.cos(deltaPhi(phi1,phi2))
        sinheta1Xsinheta2 = np.sinh(eta1)*np.sinh(eta2)
        p1dotp2 = pt1*pt2*(cosphi1phi2 + sinheta1Xsinheta2)
        invm2 = E1pE22  - pt1cosheta1_sq - pt2cosheta2_sq - 2*p1dotp2
    return np.sqrt(invm2)


def zh_helper(events):
    #deltaR = (lambda obj, obj2: ak.flatten(obj.metric_table(obj2)))
    #process_sorted = (lambda proc, dr, ind: ak.where(dr > 0.8, proc, np.nan)[ind])
    ZHCand = events.FatJetSorted[:,0]
    ak4 = events.JetGood
    bjet = events.bJetGood
    qjet = events.qJetGood
    lep = events.LeptonGood

    '''
    Comparing ZH candidate with ak4, b, q ,l || l and b
    In general, want to output arrays of variables which are sorted by dR
    '''

    # ZH and AK4
    # ak4 (in 0.8 dR to ZH cand.) btag sorted by dR, ascending with nan first
    ak4_btag_indRsort = make_dRsorted_arr(ZHCand, ak4, "lt", ak4.btagDeepFlavB)

    # ZH and b

    b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort = [make_dRsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_dRsort_vec = zip_4vec(b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort)

    b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort = [make_ptsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_ptsort_vec = zip_4vec(b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort)

    # ZH and q

    q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort = [make_dRsorted_arr(ZHCand, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass]]
    q_dRsort_vec = zip_4vec(q_pt_dRsort, q_eta_dRsort, q_phi_dRsort, q_mass_dRsort)

    q_pt_ptsort, q_eta_ptsort, q_phi_ptsort, q_mass_ptsort = [make_ptsorted_arr(ZHCand, qjet, "gt", proc) for proc in [qjet.pt, qjet.eta, qjet.phi, qjet.mass]]
    q_ptsort_vec = zip_4vec(q_pt_ptsort, q_eta_ptsort, q_phi_ptsort, q_mass_ptsort)

    # Combinations

    ZH_l_dr = deltaR(ZHCand, lep)

    ZH_l_invM = (ZHCand + lep).mass #might need this to be specifically pNet mass

    # Leading-Subleading (pt) b and q combinations

    ind_leading_b = ak.argmax(ak.nan_to_num(b_pt_ptsort, nan=-1), axis = 1, keepdims=True)

    leading_b = b_ptsort_vec[ind_leading_b]

    ind_subLeading_b = np.where(ind_leading_b > 0, -1, ind_leading_b)

    subLeading_b = b_ptsort_vec[ind_subLeading_b]

    b_b_dr = deltaR(leading_b, subLeading_b)

    ind_leading_q = ak.argmax(ak.nan_to_num(q_pt_ptsort, nan=-1), axis = 1, keepdims=True)

    leading_q = q_ptsort_vec[ind_leading_q]

    ind_subLeading_q = np.where(ind_leading_q > 0, -1, ind_leading_q)

    subLeading_q = q_ptsort_vec[ind_subLeading_q]

    q_q_dr = deltaR(leading_q, subLeading_q)

    # Nearest and second nearest b to l
    # Any NaN value is a b within the dR cone of the ZH candidate

    b_pt_dRsort_l, b_eta_dRsort_l, b_phi_dRsort_l, b_mass_dRsort_l = [make_dRsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_dRsort_l_vec = zip_4vec(b_pt_dRsort_l, b_eta_dRsort_l, b_phi_dRsort_l, b_mass_dRsort_l)

    ind_near_b = ak.argmax(ak.nan_to_num(b_pt_dRsort_l, nan=-1), axis = 1, keepdims=True)

    near_b = b_dRsort_l_vec[ind_near_b]

    ind_subNear_b = np.where(ind_near_b > 0, -1, ind_near_b)

    subNear_b = b_dRsort_l_vec[ind_subNear_b]

    l_b1_invM = (lep + near_b).mass
    l_b1_dr = deltaR(lep, near_b)

    l_b2_invM = (lep + subNear_b).mass
    l_b2_dr = deltaR(lep, subNear_b)

    # Combinations involving ZH and leading(subleading) b and q

    events["outZH_b1_pt"] = leading_b.pt


