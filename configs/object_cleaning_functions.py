import awkward as ak
import numpy as np

'''
Custom function to remove isolation bit from WPBitmap
'''

def calc_cutBasedNoIso(events):
    bmap = events.Electron.vidNestedWPBitmap
    bmapCount = ak.num(bmap)
    flatbit = ak.flatten(bmap)
    flatbit = np.array(flatbit)
    noisoStorage = (flatbit >> 0) & 7
    for i in [3, 6, 9, 12, 15, 18, 24, 27]:
        noisoStorage = np.minimum(noisoStorage, flatbit >> i & 7)
    cutbasednoiso = ak.unflatten(noisoStorage, bmapCount)
    return cutbasednoiso

'''
Customized lepton selection based on pocket-coffea default
'''

def lep_sel(events, lepton_flavour, params):

    leptons = events[lepton_flavour]
    cuts = params.skim_params[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        passes_sip3d = leptons.sip3d < cuts['sip3d']
        passes_cutBasedNoIso = calc_cutBasedNoIso(events) >= cuts['cutBasedNoIso']
        passes_iso = leptons.miniPFRelIso_all < cuts['iso']

        good_leptons = passes_eta & passes_pt & passes_sip3d & passes_cutBasedNoIso & passes_iso

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_iso = leptons.miniPFRelIso_all < cuts["iso"]
        passes_id = leptons.mediumId >= cuts['id']
        passes_sip3d = leptons.sip3d < cuts['sip3d']

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id & passes_sip3d

    return leptons[good_leptons]

'''
Soft lepton selecton
'''

def soft_lep_sel(events, lepton_flavour, params):

    leptons = events[lepton_flavour]
    cuts = params.skim_params[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt <= cuts["pt"]

    if lepton_flavour == "Electron":
        passes_cutBasedNoIso = calc_cutBasedNoIso(events) >= cuts['softCutBasedNoIso']

        good_leptons = passes_eta & passes_pt & passes_cutBasedNoIso

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_id = (leptons.mediumId == cuts['softId']) | (leptons.softId == cuts['softId'])

        good_leptons = passes_eta & passes_pt & passes_id

    return leptons[good_leptons]

'''
FatJet
'''

def fatjet_sel(events, params, leptons_collection=""):

    fatjets = events['FatJet']
    cuts = params.object_preselection['FatJet']

    passes_eta = abs(fatjets.eta) < cuts['eta']
    passes_pt = fatjets.pt > cuts['pt']
    passes_id = fatjets.jetId >= cuts['jetId']
    passes_mass = (fatjets.particleNet_mass >= cuts['mass']['low']) & (fatjets.particleNet_mass <= cuts['mass']['high'])
    if leptons_collection != "":
        dR_jets_lep = fatjets.metric_table(events[leptons_collection])
        mask_lepton_cleaning = ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1

    good_fatjets = passes_eta & passes_pt & passes_id & passes_mass

    return fatjets[good_fatjets]

'''
b-tagging
'''

def bjet_sel(Jet, params):
    return Jet[Jet.btagDeepFlavB > params.object_preselection['btag']['wp']]
