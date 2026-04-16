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
            (ak.num(events.Muon) + ak.num(events.Electron) >= 1) &
            (ak.num(events.Jet) >= 5) &
            (ak.num(events.FatJet) >= 1) &
            (events.PuppiMET.pt >= 20) )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

precut = Cut(
        name = "preSkimCut",
        params = {},
        function = precut,
)

def network_cut(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.n_ak4jets >= 5) &
            (events.n_b_outZH == 2) &
            #(events.ZH_bbvLscore >= 0.9870) &
            (events.ZH_M >= 50) &
            (events.ZH_M <= 200)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

network_cut = Cut(
        name = "network_cut",
        params = {},
        function = network_cut,
)

'''
Baseline event selection
'''

def event_selection(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.nJetGood >= 5) &
            (events.nFatJetGood >= 1) &
            (events.PuppiMET.pt > 20) &
            (events.nMuonGood + events.nElectronGood == 1)
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

def njet_cut(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.nJetGood >= 5)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

njet_cut = Cut(
        name = "njet_cut",
        params = {},
        function = njet_cut,
)

def fatjet_cut(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = ((events.nFatJetGood >= 1))
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

fatjet_cut = Cut(
        name = "fatjet_cut",
        params = {},
        function = fatjet_cut,
)

def met_cut(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.PuppiMET.pt > 20)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

met_cut = Cut(
        name = "met_cut",
        params = {},
        function = met_cut,
)

def lep_cut(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (events.nMuonGood + events.nElectronGood == 1)
            #(events.nSoftElectronGood < 2) &
            #(events.nSoftMuonGood < 2)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

lep_cut = Cut(
        name = "lep_cut",
        params = {},
        function = lep_cut,
)


def event_variation_selection(events, params, year, sample, **kwargs):
    #Event selection mask
    mask = (
            (ak.num(events.JetGood[params["jet_pt_type"]] >= 30) >= 5) &
            (events.nFatJetGood >= 1) &
            (events.MET.pt > 20) &
            (events.nMuonGood + events.nElectronGood == 1)
            #(events.nSoftElectronGood < 2) &
            #(events.nSoftMuonGood < 2)
            )
    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

# Selection for ttbar background categorization
def get_variation_selection(jet_pt_type, name=None):
    if name == None:
        name = f"jec_{jet_pt_type}"
    return Cut(name=name, params={"jet_pt_type" : jet_pt_type}, function=event_variation_selection)


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

###########################
## Implementation of JetVetoMaps 
## (Recommended for Run 3, see https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps)
## This is mostly a default pocketcoffea function, but I am customizing to use the correct jet collection aka JetGood
def get_nano_version(events, params, year):
    '''Helper function to get the nano version from the events metadata or from the default parameters.'''
    if "nano_version" in events.metadata:
        nano_version = events.metadata["nano_version"]
    else:
        if events.metadata.get("isMC", False):
            # Try to extract from the sample name
            if "NanoAODv12" in events.metadata["filename"]:
                nano_version = 12
            elif "NanoAODv15" in events.metadata["filename"]:
                nano_version = 15
            else:
                # For MC if it's not defined we take the default nano version
                #nano_version = params["default_nano_version"][year]
                nano_version = 9
        else:
            # For data if it's not defined we take the default nano version
            #nano_version = params["default_nano_version"][year]
            nano_version = 9

    return nano_version

def compute_jetId(events, jet_type, params, year):
    """
    Add (or recompute) jet ID to the jets object based on the NanoAOD version.
    Inspired by https://gitlab.cern.ch/cms-analysis/general/HiggsDNA/-/blob/master/higgs_dna/tools/jetID.py
    """
    jets = events[jet_type]
    # Get the nano version from events metadata or from default parameters
    nano_version = get_nano_version(events, params, year)
    abs_eta = abs(jets.eta)
    print(nano_version)
    # Return the existing jetId for NanoAOD versions below 12
    if nano_version < 12:
        return jets.jetId

    # For NanoAOD version 12 and above, we recompute the jet ID criteria
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV
    elif nano_version >= 12:
        if jet_type == "JetGood":
            # Default tight
            passJetIdTight = ak.where(
                abs_eta <= 2.7,
                (jets.jetId & (1 << 1)) > 0,  # Tight criteria for abs_eta <= 2.7
                ak.where(
                    (abs_eta > 2.7) & (abs_eta <= 3.0),
                    ((jets.jetId & (1 << 1)) > 0) & (jets.neHEF < 0.99),  # Tight criteria for 2.7 < abs_eta <= 3.0
                    ((jets.jetId & (1 << 1)) > 0) & (jets.neEmEF < 0.4)  # Tight criteria for 3.0 < abs_eta
                )
            )
            # Default tight lepton veto
            passJetIdTightLepVeto = ak.where(
                abs_eta <= 2.7,
                passJetIdTight & (jets.muEF < 0.8) & (jets.chEmEF < 0.8),  # add lepton veto for abs_eta <= 2.7
                passJetIdTight  # No lepton veto for 2.7 < abs_eta
            )
            return (passJetIdTight * (1 << 1)) | (passJetIdTightLepVeto * (1 << 2))

        elif jet_type == "FatJet":
            # For nanoAOD v12, only using the original branch in the tracker acceptance
            passJetIdTight = ak.where(
                abs_eta <= 2.7,
                (jets.jetId & (1 << 1)) > 0,  # Tight criteria for abs_eta <= 2.7
                ak.zeros_like(jets.jetId) 
                )
            # Not tight lepveto for FatJet
            return (passJetIdTight * (1 << 1)) 
        else:
            raise ValueError(f"Jet type {jet_type} not recognized for JetID")


def get_JetVetoMap_custom(name="JetVetoMaps"):
    return Cut(
        name=name, params={}, function=get_JetVetoMap_Mask
    )
       
def get_JetVetoMap_Mask(events, params, year, processor_params, sample, isMC, **kwargs):
    # Import here to prevent circular import configurator -> cuts -> cut_functions -> jets -> utils -> configurator
    # For nanoV15 no jetId key in Nano anymore. 
    # For nanoV12 (i.e. 22/23), jet Id is also buggy, should therefore be rederived
    # in the following, if nano_version not explicitly specified in params, v9 is assumed for Run2UL, v12 for 22/23 and v15 for 2024
    jets = ak.with_field(events["JetGood"], compute_jetId(events, "JetGood", processor_params, year), "jetId_corrected")
    mask_for_VetoMap = (
        (jets["jetId_corrected"]>=6) # Must fulfill tightLepVeto
        & (abs(jets.eta) < 5.19) # Must be within HCal acceptance
        & (jets.pt > 15.) # Minimum pT
        & ((jets["neEmEF"]+jets["chEmEF"])<0.9) # Energy fraction not dominated by ECal
    )
    jets = jets[mask_for_VetoMap]
    cset = correctionlib.CorrectionSet.from_file(
        processor_params.jet_scale_factors.vetomaps[year]["file"]
    )
    corr = cset[processor_params.jet_scale_factors.vetomaps[year]["name"]]
    etaFlat, phiFlat, etaCounts = ak.flatten(jets.eta), ak.flatten(jets.phi), ak.num(jets.eta)
    phiFlat = np.clip(phiFlat, -3.14159, 3.14159) # Needed since no overflow included in phi binning
    weight = ak.unflatten(
        corr.evaluate("jetvetomap", etaFlat, phiFlat),
        counts=etaCounts,
    )
    eventMask = ak.sum(weight, axis=-1)==0 # if at least one jet is vetoed, reject it event
    return ak.where(ak.is_none(eventMask), False, eventMask)

'''
Helper function to sort jets by highest Xbb score
'''

def sortbyscore(coll, score):
    return coll[ak.argsort(coll[score], axis=1, ascending=False)]
