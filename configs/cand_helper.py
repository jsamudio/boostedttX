import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector
# We do this for various objects so why not make it a single function

deltaR = (lambda obj1_, obj2_: ak.flatten(obj1_.metric_table(obj2_)))

def make_dRsorted_arr(obj1, obj2, comparison, proc):
    if (comparison == "gt"): #greaterthan
        nanArgsort = (lambda dr_: ak.argsort(ak.where(dr_ > 0.8, dr_, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ > 0.8, proc_, np.nan)[ind_])
    else:
        nanArgsort = (lambda dr_: ak.argsort(ak.where(dr_ <= 0.8, dr_, np.nan), axis=1))
        process_sorted = (lambda proc_, dr_, ind_: ak.where(dr_ <= 0.8, proc_, np.nan)[ind_])

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


def zh_helper(events):
    #deltaR = (lambda obj, obj2: ak.flatten(obj.metric_table(obj2)))
    #process_sorted = (lambda proc, dr, ind: ak.where(dr > 0.8, proc, np.nan)[ind])
    ZHCand = events.FatJetSorted[:,0]
    ak4 = events.JetGood
    bjet = events.bJetGood
    '''
    Comparing ZH candidate with ak4, b, q ,l || l and b
    In general, want to output arrays of variables which are sorted by dR
    '''

    # ZH and AK4
    # ak4 (in 0.8 dR to ZH cand.) btag sorted by dR, ascending with nan first
    ak4_btag_indRsort = make_dRsorted_arr(ZHCand, ak4, "lt", ak4.btagDeepFlavB)

    #ZH and b

    b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort = [make_dRsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_dRsort_vec = zip_4vec(b_pt_dRsort, b_eta_dRsort, b_phi_dRsort, b_mass_dRsort)

    b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort = [make_ptsorted_arr(ZHCand, bjet, "gt", proc) for proc in [bjet.pt, bjet.eta, bjet.phi, bjet.mass]]
    b_ptsort_vec = zip_4vec(b_pt_ptsort, b_eta_ptsort, b_phi_ptsort, b_mass_ptsort)

    print(b_pt_ptsort)
