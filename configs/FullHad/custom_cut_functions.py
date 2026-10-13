import awkward as ak
import numpy as np
from pocket_coffea.lib.cut_definition import Cut
from collections.abc import Iterable

'''
Precut to slim processing
'''

def precut(events, params, year, sample, **kwargs):
    mask = (
            (ak.num(events.Jet) >= 6) &
            (ak.num(events.FatJet) >= 1))
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
    mask = (
            (events.nJetGood >= 6) &
            (events.nFatJetGood >= 1) &
            (ak.sum(events.Jet.pt) > 500) &
            (events.nMuonGood + events.nElectronGood == 0)
            )
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
    mask = (
            (events.nJetGood >= 5) &
            (events.nFatJetGood >= 1) &
            (events.MET.pt > 20) &
            (events.nMuonGood + events.nElectronGood == 1)
            )
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
─────────────────────────────────────────────────────────────
Affected detector region:
    -3.0 < eta < -1.3  AND  -1.57 < phi < -0.87

Two HCAL modules (HEM15 and HEM16) failed during Run2018.
Affected runs: run >= 319077 (part of Run2018C + all of Run2018D)
Affected lumi fraction: ~64.9% of total 2018 lumi (38.77 / 59.74 fb-1)

Data: veto events with a jet/fatjet in HEM region for run >= 319077
MC:   veto the same fraction of events randomly, using a deterministic
      per-event hash so the same event is always treated the same way
      regardless of chunk boundaries.
Normalisation: unchanged — lumi weight always uses full 2018 lumi.
─────────────────────────────────────────────────────────────
'''

# Lumi fraction of HEM-affected runs in 2018
HEM_LUMI_FRACTION = 0.649

def HEM_veto_func(events, params, year, sample, **kwargs):

    # Only relevant for 2018
    if year != "2018":
        return ak.ones_like(events.event, dtype=bool)

    # HEM affected detector region
    hem_eta_min, hem_eta_max = -3.0, -1.3
    hem_phi_min, hem_phi_max = -1.57, -0.87

    # Check AK4 jets in HEM region (pt > 15 GeV)
    jets = events.Jet
    jet_in_hem = (
        (jets.eta > hem_eta_min) & (jets.eta < hem_eta_max) &
        (jets.phi > hem_phi_min) & (jets.phi < hem_phi_max) &
        (jets.pt > 15)
    )

    # Check AK8 fatjets in HEM region (pt > 15 GeV)
    fatjets = events.FatJet
    fatjet_in_hem = (
        (fatjets.eta > hem_eta_min) & (fatjets.eta < hem_eta_max) &
        (fatjets.phi > hem_phi_min) & (fatjets.phi < hem_phi_max) &
        (fatjets.pt > 15)
    )

    # Event has at least one object in HEM region
    has_hem_object = (
        ak.any(jet_in_hem, axis=1) |
        ak.any(fatjet_in_hem, axis=1)
    )

    isMC = str(events.metadata.get("isMC", "False")).lower() in ("true", "1")

    if isMC:
        # Deterministic per-event random decision using Knuth multiplicative hash
        # Same event number → same veto decision regardless of chunk
        event_hash = (ak.to_numpy(events.event) * 2654435761) & 0xFFFFFFFF
        in_affected_period = ak.Array(event_hash % 1000 < int(HEM_LUMI_FRACTION * 1000))
    else:
        # Data: only veto events in runs after HEM failure
        in_affected_period = events.run >= 319077

    # Veto: event is in affected period AND has an object in HEM region
    veto = ~(in_affected_period & has_hem_object)

    return ak.where(ak.is_none(veto), True, veto)

HEM_veto = Cut(
    name="HEM_veto",
    params={},
    function=HEM_veto_func,
)
