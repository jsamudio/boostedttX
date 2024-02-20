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
            (events.MET.pt > 20) &
            (events.nElectronGood + events.nMuonGood == 1) )
    # Pad None vlaues with False
    return ak.where(ak.is_none(mask), False, mask)

event_selection = Cut(
        name = "eventSelection",
        params = {},
        function = event_selection,
)

def btag_mask(events, params, year, sample, **kwargs):
    mask = ((events.nbJetGood >= 2))
    return ak.where(ak.is_none(mask), False, mask)

btag_mask = Cut(
        name = "btagMask",
        params = {},
        function = btag_mask,
)
