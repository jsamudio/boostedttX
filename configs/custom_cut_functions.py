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
        name = "eventSelection",
        params = {},
        function = event_selection,
)

def calc_cutBasedNoIso(events, params):
    bmap = events['Electron']['Electron_vidNestedWPBitmap']
    bmapCount = bmap.counts
    flatbit = bmap.flatten()
    flatbit = np.array(flatbit)
    noiso_storage = (flatbit >> 0) & 7
    for i in [3, 6, 9, 12, 15, 18, 24, 27]:
        noiso_storage = np.minimum(noiso_storage, flatbit >> i & 7)
    cutbasednoiso = ak.unflatten(noiso_storage, bmapcount)
    return cutbasednoiso

def lep_sel(events, lepton_flavour, params):

    leptons = events[lepton_flavour]
    cuts = params.object_preselection[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        passes_sip3d = leptons.sip3d < cuts['sip3d']
        passes_cutbasedNoIso = leptons.cutBasedNoIso >= cuts['cutBasedNoIso']
        passes_iso = leptons.miniPFRelIso_all < cuts['iso']

        good_leptons = passes_eta & passes_pt & passes_sip3d & passes_cutBasedNoIso & passes_iso

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_iso = leptons.miniPFRelIso_all < cuts["iso"]
        passes_id = leptons.mediumId >= cuts['id']
        passes_sip3d = leptons.sip3d < cuts['sip3d']

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id & passes_sip3d

    return leptons[good_leptons]
