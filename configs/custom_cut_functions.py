import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

def precut(events, params, year, sample, **kwargs):
    #Precut Mask
    mask = (
            (events['nMuon'] + events['nElectron'] >= 1) &
            (events['nJet'] >= 5) &
            (events['nFatJet'] >= 1) &
            (events['MET_pt'] >= 20) )
    # Pad None values with False
    return ak.where(ak.is_one(mask), False, mask)

precut = Cut(
        name = "preSkimCut",
        params = {},
        function = precut,
)

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
        name = "eventSelection"
        params = {},
        function = event_selection,
)

###FIXME customize what is a good lepton since we use some additional variables.
def lepton_selection(events, lepton_flavour, params):

    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        # Requirements on SuperCluster eta, isolation and id
        etaSC = abs(leptons.deltaEtaSC + leptons.eta)
        passes_SC = np.invert((etaSC >= 1.4442) & (etaSC <= 1.5660))
        passes_iso = leptons.pfRelIso03_all < cuts["iso"]
        passes_id = leptons[cuts['id']] == True

        good_leptons = passes_eta & passes_pt & passes_SC & passes_iso & passes_id

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_iso = leptons.pfRelIso04_all < cuts["iso"]
        passes_id = leptons[cuts['id']] == True

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id

    return leptons[good_leptons]
