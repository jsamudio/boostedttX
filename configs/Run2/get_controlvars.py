import math
import awkward as ak
import numpy as np
from coffea.nanoevents.methods import vector

FILL = -1.0   # value for events that have fewer objects than requested


def sortPt(obj):
    """Objects sorted by pT, descending."""
    return obj[ak.argsort(obj.pt, ascending=False)]


def nth(objs, n, field, fill=FILL):
    """`field` of the n-th object per event (0-based); `fill` where the event has fewer objects.
    Needed since the b-tag requirement moved from the preselection to the btag_mask category:
    events with fewer jets than expected now reach this code."""
    padded = ak.pad_none(objs, n + 1, clip=True)
    return ak.fill_none(padded[:, n][field], fill)


def get_controlvars(events):
    ak4 = sortPt(events.JetGood)
    ak8 = sortPt(events.FatJetGood)
    bjet = sortPt(events.bJetGood)
    qjet = events.qJetGood
    lep = events.LeptonGood
    ele = events.ElectronGood
    muon = events.MuonGood
    met = events.MET          # type-1 corrected PF MET (Run 2), as in the selection and SPANet inputs

    # leading is 1 and subleading is 2

    events['nPV'] = events.PV.npvs
    events['nPVGood'] = events.PV.npvsGood

    events['MET_pt'] = met.pt
    events['MET_phi'] = met.phi

    events['lep_pt'] = ak.firsts(lep.pt)
    events['ele_pt'] = ak.firsts(ele.pt)
    events['muon_pt'] = ak.firsts(muon.pt)
    events['lep_eta'] = ak.firsts(lep.eta)
    events['ele_eta'] = ak.firsts(ele.eta)
    events['muon_eta'] = ak.firsts(muon.eta)

    # Relative electron energy shifts of the leading selected electron (EGM scale / smearing check,
    # ele_energy_check.py); 0 without an electron or when the NanoAOD has no such field
    ele0 = ak.firsts(ele)
    energy = ele0.pt * np.cosh(ele0.eta)
    for var in ('dEscaleUp', 'dEscaleDown', 'dEsigmaUp', 'dEsigmaDown'):
        rel = ele0[var] / energy if var in ele.fields else ak.zeros_like(ele0.pt)
        events[f'ele_{var}'] = ak.fill_none(rel, 0.0)

    events['n_ak4'] = ak.count(ak4.pt, axis=1)
    events['n_bjet'] = ak.count(bjet.pt, axis=1)
    events['n_ak8'] = ak.count(ak8.pt, axis=1)

    events['jet1_pt'] = nth(ak4, 0, 'pt')
    events['jet2_pt'] = nth(ak4, 1, 'pt')
    events['bjet1_pt'] = nth(bjet, 0, 'pt')
    # events['bjet2_pt'] = nth(bjet, 1, 'pt')        # safe to enable now

    events['jet1_eta'] = nth(ak4, 0, 'eta')
    events['jet2_eta'] = nth(ak4, 1, 'eta')
    events['bjet1_eta'] = nth(bjet, 0, 'eta')
    # events['bjet2_eta'] = nth(bjet, 1, 'eta')      # safe to enable now

    events['jet1_btag'] = nth(ak4, 0, 'btagB')
    events['jet2_btag'] = nth(ak4, 1, 'btagB')
    events['bjet1_btag'] = nth(bjet, 0, 'btagB')
    # events['bjet2_btag'] = nth(bjet, 1, 'btagB')   # safe to enable now (was the FIXME)

    events['fatjet1_pt'] = nth(ak8, 0, 'pt')
    events['fatjet1_eta'] = nth(ak8, 0, 'eta')
    events['fatjet1_mass'] = nth(ak8, 0, 'mass')