import awkward as ak
import numpy as np

#customized for boosted ZHbb

def lepton_selection(events, lepton_flavour, params):

    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]
    # pT and eta cuts
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        passed_sip3d = leptons.sip3d < cuts["sip3d"]
        passed_cutBasedNoIso = leptons.cutBasedNoIso < cuts["cutBasedNoIso"]
        passed_iso = leptons.miniPFRelIso_all < cuts["miniPFRelIso_all"]

    good_leptons = passes_eta & passes_pt & passes_sip3d & passes_cutBasedNoIso & passed_iso
