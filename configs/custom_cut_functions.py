import awkward as ak
import numpy as np
from pocket_coffea.lib.cut_definition import Cut

'''
Precut to slim processing
'''

def precut(events, params, year, sample, **kwargs):
    #Precut Mask
    mask = (
            (ak.num(events.Muon) + ak.num(events.Electron) >= 1) &
            (ak.num(events.Jet) >= 5) &
            (ak.num(events.FatJet) >= 1) &
            (events.MET.pt >= 20) )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

precut = Cut(
        name = "preSkimCut",
        params = {},
        function = precut,
)

'''
Baseline event selection
'''

def event_selection(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.nJetGood >= 5) &
            (events.nFatJetGood >= 1) &
            (events['MET_pt'] > 20) &
            (events.nGoodElectron + events.nGoodMuon == 1) )
    # Pad None vlaues with False
    return ak.where(ak.is_one(mask), False, mask)

event_selection = Cut(
        name = "eventSelection",
        params = {},
        function = event_selection,
)

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
