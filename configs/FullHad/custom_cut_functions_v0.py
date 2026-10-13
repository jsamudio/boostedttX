import awkward as ak
import numpy as np
from pocket_coffea.lib.cut_definition import Cut
from collections.abc import Iterable

'''
Precut to slim processing
'''

def precut(events, params, year, sample, **kwargs):
    #Precut Mask
    mask = (
            #(ak.num(events.Muon) + ak.num(events.Electron) == 0) &
            (ak.num(events.Jet) >= 6) &
            (ak.num(events.FatJet) >= 1))
            #(events.MET.pt >= 20) )
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
            (events.nJetGood >= 6) &
            (events.nFatJetGood >= 1) &
            (ak.sum(events.Jet.pt) > 500) &
            #(events.MET.pt > 20) &
            (events.nMuonGood + events.nElectronGood == 0)
            #(events.nSoftElectronGood < 2) &
            #(events.nSoftMuonGood < 2)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

event_selection = Cut(
        name = "eventSelection",
        params = {},
        function = event_selection,
)

'''
Leptons study event selection
'''

def lepstudy_event_selection(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.nJetGood >= 5) &
            (events.nFatJetGood >= 1) &
            (events.MET.pt > 20) &
            (events.nMuonGood + events.nElectronGood == 1)
            #(events.nSoftElectronGood < 2) &
            #(events.nSoftMuonGood < 2)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

lepstudy_event_selection = Cut(
        name = "lepstudy_eventSelection",
        params = {},
        function = lepstudy_event_selection,
)

'''
btag mask
'''

def btag_mask(events, params, year, sample, **kwargs):
    mask = ((events.nbJetGood >= 2))
    return ak.where(ak.is_none(mask), False, mask)

btag_mask = Cut(
        name = "btagMask",
        params = {},
        function = btag_mask,
)

'''
btag mask
'''

#def eleMVAIsoWPL_mask(events, params, year, sample, **kwargs):
#    mask = ((events.nMuonGood + events.nElectronGood_mvaIsoWPL == 1))
#    return ak.where(ak.is_none(mask), False, mask)

#btag_mask = Cut(
#        name = "ele_isoWPL",
#        params = {},
#        function = btag_mask,
#)

'''
SFOS J/Psi veto
'''
def vetoMu(events, params, year, sample, **kwargs):
    OS = events.mu_softmu.charge == 0

    mask = (
        OS &
        (events.mu_softmu.mass < 12) &
        (events.mu_softmu.mass > 0)
        )
    return ak.where(ak.is_none(mask), False, ~mask)

vetoMu = Cut(
        name= "vetoMu",
        params = {},
        function = vetoMu,
)

def vetoE(events, params, year, sample, **kwargs):
    OS = events.e_softe.charge == 0

    mask = (
        OS &
        (events.e_softe.mass < 12) &
        (events.e_softe.mass > 0)
        )
    return ak.where(ak.is_none(mask), False, ~mask)

vetoE = Cut(
        name= "vetoE",
        params = {},
        function = vetoE,
)

'''
cut function for tt
'''

def eq_genTtbarId_100(events, params, year, sample, **kwargs):
    """
    This function returns a mask for events where genTtbarId % 100 == params["genTtbarId"]
    or the logic OR of masks in the case in which params["genTtbarId"] is an iterable.
    The dictionary for genTtbarId % 100 is the following
    (taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/GenHFHadronMatcher, visited on 23.11.2023):
    0  : "tt+LF",
    41 : "tt+c",
    42 : "tt+2c",
    43 : "tt+cc",
    44 : "tt+c2c",
    45 : "tt+2c2c",
    46 : "tt+C",
    51 : "tt+b",
    52 : "tt+2b",
    53 : "tt+bb",
    54 : "tt+b2b",
    55 : "tt+2b2b",
    56 : "tt+B",
    """
    allowed_ids = [0, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56]
    if type(params["genTtbarId"]) == int:
        if params["genTtbarId"] in allowed_ids:
            return events.genTtbarId % 100 == params["genTtbarId"]
        else:
            raise Exception(f"The cut on genTtbarId % 100 must be an integer between 0 and 56.\nPossible choices:{allowed_ids}")
    elif isinstance(params["genTtbarId"], Iterable):
        mask = ak.zeros_like(events.event, dtype=bool)
        for _id in params["genTtbarId"]:
            if _id in allowed_ids:
                mask = mask | (events.genTtbarId % 100 == _id)
            else:
                raise Exception(f"The cut on genTtbarId % 100 must be an integer between 0 and 56.\nPossible choices:{allowed_ids}")
        return mask
    else:
        raise Exception(f'params["genTtbarId"] must be an integer or an iterable of integers between 0 and 56.\nPossible choices:{allowed_ids}')

# Selection for ttbar background categorization
def get_genTtbarId_100_eq(genTtbarId, name=None):
    if name == None:
        if type(genTtbarId) == int:
            name = f"genTtbarId_100_eq_{genTtbarId}"
        if isinstance(genTtbarId, Iterable):
            name = f"genTtbarId_100_eq_" + "_".join([str(s) for s in genTtbarId])
    return Cut(name=name, params={"genTtbarId" : genTtbarId}, function=eq_genTtbarId_100)

'''
GenMatch Cut
'''

def genMatchZHbb(events, params, year, sample, **kwargs):
    mask = (events['matchedGen_ZHbb_bb'])

    return ak.where(ak.is_none(mask), False, mask)

def non_genMatchZHbb(events, params, year, sample, **kwargs):
    mask = (events['matchedGen_ZHbb_bb'])

    return ak.where(ak.is_none(mask), False, ~mask)

genMatch = Cut(
        name = "genMatch",
        params = {},
        function = genMatchZHbb,
)

non_genMatch = Cut(
        name = "genMatch",
        params = {},
        function = non_genMatchZHbb,
)

'''
Helper function to sort jets by highest Xbb score
'''

def sortbyscore(coll, score):
    return coll[ak.argsort(coll[score], axis=1, ascending=False)]
'''
HEM15/16 veto for 2018
Affected region: -3.0 < eta < -1.3 and -1.57 < phi < -0.87
Affects Run2018C (partial) and Run2018D (run >= 319077)
MC: apply randomly to 64.9% of events (lumi fraction of affected runs)
'''

def HEM_veto(events, params, year, sample, **kwargs):
    # Only apply for 2018
    if year != "2018":
        return ak.ones_like(events.event, dtype=bool)
    
    # HEM affected detector region
    hem_eta_min, hem_eta_max = -3.0, -1.3
    hem_phi_min, hem_phi_max = -1.57, -0.87

    # Check if any jet or FatJet falls in HEM region
    jets = events.Jet
    fatjets = events.FatJet

    jet_in_hem = (
        (jets.eta > hem_eta_min) & (jets.eta < hem_eta_max) &
        (jets.phi > hem_phi_min) & (jets.phi < hem_phi_max) &
        (jets.pt > 15)
    )
    fatjet_in_hem = (
        (fatjets.eta > hem_eta_min) & (fatjets.eta < hem_eta_max) &
        (fatjets.phi > hem_phi_min) & (fatjets.phi < hem_phi_max) &
        (fatjets.pt > 15)
    )

    has_hem_object = (
        ak.any(jet_in_hem, axis=1) |
        ak.any(fatjet_in_hem, axis=1)
    )

    isMC = events.metadata.get("isMC", "False") == "True"

    if isMC:
        # Apply to 64.9% of MC events (lumi fraction: affected/total 2018)
        # Use event number as seed for reproducibility
        rng = np.random.default_rng(seed=42)
        random_mask = ak.Array(rng.random(len(events)) < 0.649)
        run_mask = random_mask
    else:
        # For data: only veto events in runs >= 319077
        run_mask = events.run >= 319077

    # Veto events that have a HEM object AND are in the affected run range
    veto = ~(run_mask & has_hem_object)

    return ak.where(ak.is_none(veto), True, veto)

HEM_veto = Cut(
    name="HEM_veto",
    params={},
    function=HEM_veto,
)
