"""Object selection of the Run 3 analysis (cuts in params/object_preselection.yaml)."""
import awkward as ak
import numpy as np
from pocket_coffea.lib.jets import compute_jetId


def lep_sel(events, lepton_flavour, params):
    """Tight leptons: electrons MVA noIso WP90, muons medium ID; mini-isolation and sip3d for both."""
    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]
    good = (abs(leptons.eta) < cuts["eta"]) & (leptons.pt > cuts["pt"]) & (leptons.sip3d < cuts["sip3d"])
    good = good & (leptons.miniPFRelIso_all < cuts["iso"])
    if lepton_flavour == "Electron":
        good = good & (leptons[cuts["id"]] == True)
    elif lepton_flavour == "Muon":
        good = good & (leptons.mediumId >= cuts["id"])
    return leptons[good]


def soft_lep_sel(events, lepton_flavour, params):
    """Soft leptons (below the tight pT threshold) for the J/psi-like OS pair veto."""
    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]
    good = (abs(leptons.eta) < cuts["eta"]) & (leptons.pt <= cuts["pt"])
    if lepton_flavour == "Electron":
        good = good & (leptons[cuts["id"]] == True)
    elif lepton_flavour == "Muon":
        good = good & ((leptons.mediumId == cuts["softId"]) | (leptons.softId == cuts["softId"]))
    return leptons[good]


def lep_softlep_combo(lepton, softLepton, transverse=False):
    """Sum of the leading tight/soft lepton pair (or zeros if fewer than two leptons)."""
    fields = {"pt": 0., "eta": 0., "phi": 0., "mass": 0., "charge": 0.}
    leptons = ak.pad_none(ak.with_name(ak.concatenate([lepton[:, 0:2], softLepton[:, 0:2]], axis=1),
                                       "PtEtaPhiMCandidate"), 2)
    nlep = ak.num(leptons[~ak.is_none(leptons, axis=1)])
    ll = leptons[:, 0] + leptons[:, 1]
    for var in fields.keys():
        fields[var] = ak.where((nlep == 2), getattr(ll, var), fields[var])
    fields["deltaR"] = ak.where((nlep == 2), leptons[:, 0].delta_r(leptons[:, 1]), -1)
    if transverse:
        fields["eta"] = ak.zeros_like(fields["pt"])
    return ak.zip(fields, with_name="PtEtaPhiMCandidate")


AK4_TAGGERS = {
    "UParT": ("btagUParTAK4B", "btagUParTAK4CvL", "btagUParTAK4CvB"),
    "PNet": ("btagPNetB", "btagPNetCvL", "btagPNetCvB"),
    "DeepFlav": ("btagDeepFlavB", "btagDeepFlavCvL", "btagDeepFlavCvB"),
}


def jet_selection_custom(events, jet_type, params, year, leptons_collection="", jet_tagger=""):
    """AK4 (`Jet`) or AK8 (`FatJet`) selection: pT, eta, jet ID (recomputed for NanoAODv15),
    lepton cleaning; AK8 also the mass window on FatJet.mass (the GloParT regressed mass, set in
    the workflow). Adds btagB/btagCvL/btagCvB (AK4) or btagBB/btagCC = GloParT X/(X+QCD) (AK8)."""
    jets = events[jet_type]
    cuts = params.object_preselection[jet_type]
    jets["jetId_corrected"] = compute_jetId(events, jet_type, params, year)
    mask = (jets.pt > cuts["pt"]) & (np.abs(jets.eta) < cuts["eta"]) & (jets.jetId_corrected >= cuts["jetId"])
    if leptons_collection != "":
        dR_jets_lep = jets.metric_table(events[leptons_collection])
        mask = mask & (ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1)

    if jet_type == "Jet":
        if jet_tagger:
            key = next((k for k in AK4_TAGGERS if k in jet_tagger), None)
            if key is None:
                raise NotImplementedError(f"This tagger is not implemented: {jet_tagger}")
            B, CvL, CvB = AK4_TAGGERS[key]
            jets["btagB"], jets["btagCvL"], jets["btagCvB"] = jets[B], jets[CvL], jets[CvB]

    elif jet_type == "FatJet":
        mask = mask & (jets.mass > cuts["mass_min"]) & (jets.mass < cuts["mass_max"])
        if jet_tagger:
            if "GloParT" not in jet_tagger:
                raise NotImplementedError(f"This tagger is not implemented: {jet_tagger}")
            jets["btagBB"] = jets.globalParT3_Xbb / (jets.globalParT3_Xbb + jets.globalParT3_QCD)
            jets["btagCC"] = jets.globalParT3_Xcc / (jets.globalParT3_Xcc + jets.globalParT3_QCD)

    return jets[mask], mask


def btag_wp(params, year):
    """AK4 b-tagging threshold of the working point in object_preselection.btag.wp (e.g. 'M')."""
    return params.btagging.working_point[year]["btagging_WP"][params.object_preselection.btag.wp]


def bjet_sel(Jet, params, year):
    return Jet[Jet.btagB > btag_wp(params, year)]


def qjet_sel(Jet, params, year):
    return Jet[Jet.btagB <= btag_wp(params, year)]
