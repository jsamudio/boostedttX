import awkward as ak
import numpy as np
from pocket_coffea.lib.jets import compute_jetId

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
    cuts = params.object_preselection[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        passes_sip3d = leptons.sip3d < cuts['sip3d']
        #passes_cutBasedNoIso = calc_cutBasedNoIso(events) >= cuts['cutBasedNoIso']
        passes_id = leptons.mvaNoIso_WP90 == True
        passes_iso = leptons.miniPFRelIso_all < cuts['iso']

        good_leptons = passes_eta & passes_pt & passes_sip3d & passes_id & passes_iso

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_iso = leptons.miniPFRelIso_all < cuts["iso"]
        passes_id = leptons.mediumId >= cuts['id']
        passes_sip3d = leptons.sip3d < cuts['sip3d']

        good_leptons = passes_eta & passes_pt & passes_iso & passes_id & passes_sip3d

    return leptons[good_leptons]

def lepstudy_sel(events, lepton_flavour, params, electron_id='mvaFall17V2noIso_WP90'):

    leptons = events[lepton_flavour]
    cuts = params.skim_params[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt > cuts["pt"]

    if lepton_flavour == "Electron":
        passes_sip3d = leptons.sip3d < cuts['sip3d']
        passes_id = leptons[electron_id] == True
        if 'noIso' in electron_id:
            passes_iso = leptons.miniPFRelIso_all < cuts['iso']
        else:
            passes_iso = True

        good_leptons = passes_eta & passes_pt & passes_sip3d & passes_id & passes_iso

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
    cuts = params.object_preselection[lepton_flavour]
    # Requirements on pT and eta
    passes_eta = abs(leptons.eta) < cuts["eta"]
    passes_pt = leptons.pt <= cuts["pt"]

    if lepton_flavour == "Electron":
        passes_id = leptons.mvaNoIso_WP90 == True

        good_leptons = passes_eta & passes_pt & passes_id

    elif lepton_flavour == "Muon":
        # Requirements on isolation and id
        passes_id = (leptons.mediumId == cuts['softId']) | (leptons.softId == cuts['softId'])

        good_leptons = passes_eta & passes_pt & passes_id

    return leptons[good_leptons]

'''
Get lep + SFOS lep
'''
def lep_softlep_combo(lepton, softLepton, transverse=False):
    fields = {
        "pt": 0.,
        "eta": 0.,
        "phi": 0.,
        "mass": 0.,
        "charge": 0.,
    }

    leptons = ak.pad_none(ak.with_name(ak.concatenate([ lepton[:, 0:2], softLepton[:, 0:2]], axis=1), "PtEtaPhiMCandidate"), 2)
    nlep =  ak.num(leptons[~ak.is_none(leptons, axis=1)])
    ll = leptons[:,0] + leptons[:,1]

    for var in fields.keys():
        fields[var] = ak.where(
            (nlep == 2),
            getattr(ll, var),
            fields[var]
        )

    fields["deltaR"] = ak.where(
        (nlep == 2), leptons[:,0].delta_r(leptons[:,1]), -1)

    if transverse:
        fields["eta"] = ak.zeros_like(fields["pt"])
    dileptons = ak.zip(fields, with_name="PtEtaPhiMCandidate")

    return dileptons


'''
FatJet
'''

def ak_crosscleaned(eta1,phi1,eta2,phi2,cut):
    #input 4 awkward arrays
    #jet is second
    import awkward as ak
    import math
    c1_ = np.array(ak.count(eta1, axis=-1))
    c2_ = np.array(ak.count(eta2, axis=-1))
    out_ = np.ones(len(ak.flatten(eta2)))
    args_ = ak.flatten(eta1), ak.flatten(phi1), c1_, ak.flatten(eta2), ak.flatten(phi2), c2_, cut, math.pi, out_
    def cc(e1,p1,c1,e2,p2,c2,cut,pi,o):
        iter1 = 0
        iter2 = 0
        a = c1.size
        for i in range(a):
            for ij in range(c2[i]):
                for ii in range(c1[i]):
                    deta = e1[ii+iter1] - e2[ij+iter2]
                    dphi = p1[ii+iter1] - p2[ij+iter2]
                    if (dphi > pi):
                        dphi = dphi - 2*pi
                    elif (dphi <= -pi):
                        dphi = dphi + 2*pi
                    dr = float((deta**2 + dphi**2)**(.5))
                    if dr <= cut:
                        o[ij+iter2] = 0
                        break
            iter1 += c1[i]
            iter2 += c2[i]
        return o
    out = ak.unflatten(cc(*args_), ak.count(eta2,axis=-1))
    out = ak.values_astype(out, bool)
    return out

def is_lep_cleaned(events, lepton_collection, jet_collection, cut=0.4):
    lep = events[lepton_collection]
    jets = events[jet_collection]
    lep_eta, lep_phi = lep.eta, lep.phi
    jets_eta, jets_phi = jets.eta, jets.phi
    jet_mask = ak_crosscleaned(lep_eta, lep_phi, jets_eta, jets_phi, cut)
    return jet_mask

def jet_selection_custom(events, jet_type, params, year, leptons_collection="", jet_tagger=""):
    jets = events[jet_type]
    cuts = params.object_preselection[jet_type]

    # For nanoV15 no jetId key in Nano anymore. 
    # For nanoV12 (i.e. 22/23), jet Id is also buggy, should therefore be rederived
    # in the following, if nano_version not explicitly specified in params, v9 is assumed for Run2UL, v12 for 22/23 and v15 for 2024
    jets["jetId_corrected"] = compute_jetId(events, jet_type, params, year)

    # Mask for  jets not passing the preselection
    mask_presel = (
        (jets.pt > cuts["pt"])
        & (np.abs(jets.eta) < cuts["eta"])
        & (jets.jetId_corrected >= cuts["jetId"])
    )
    # Lepton cleaning
    # Only jets that are more distant than dr to ALL leptons are tagged as good jets
    if leptons_collection != "":
        dR_jets_lep = jets.metric_table(events[leptons_collection])
        mask_lepton_cleaning = ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1
    else:
        mask_lepton_cleaning = True

    if jet_type == "Jet":
        # Selection on PUid. Only available in Run2 UL, thus we need to determine which sample we run over;
        if year in ['2016_PreVFP', '2016_PostVFP','2017','2018']:
            mask_jetpuid = (jets.puId >= params.jet_scale_factors.jet_puId[year]["working_point"][cuts["puId"]["wp"]]) | (
                jets.pt >= cuts["puId"]["maxpt"]
            )
        else:
            mask_jetpuid = True
  
        mask_good_jets = mask_presel & mask_lepton_cleaning & mask_jetpuid

        if jet_tagger != "":
            if "PNet" in jet_tagger:
                B   = "btagPNetB"
                CvL = "btagPNetCvL"
                CvB = "btagPNetCvB"
            elif "DeepFlav" in jet_tagger:
                B   = "btagDeepFlavB"
                CvL = "btagDeepFlavCvL"
                CvB = "btagDeepFlavCvB"
            elif "RobustParT" in jet_tagger:
                B   = "btagRobustParTAK4B"
                CvL = "btagRobustParTAK4CvL"
                CvB = "btagRobustParTAK4CvB"
            elif "UParT" in jet_tagger:
                B   = "btagUParTAK4B"
                CvL = "btagUParTAK4CvL"
                CvB = "btagUParTAK4CvB"
            else:
                raise NotImplementedError(f"This tagger is not implemented: {jet_tagger}")
            
            if B not in jets.fields or CvL not in jets.fields or CvB not in jets.fields:
                raise NotImplementedError(f"{B}, {CvL}, and/or {CvB} are not available in the input.")

            jets["btagB"] = jets[B]
            jets["btagCvL"] = jets[CvL]
            jets["btagCvB"] = jets[CvB]

    elif jet_type == "FatJet":
        # Apply the msd and preselection cuts
        regressed_mass = events.FatJet.globalParT3_massCorrX2p * events.FatJet.mass * (1 - events.FatJet.rawFactor)
        mask_msd = (regressed_mass > cuts["mass_min"]) & (regressed_mass < cuts["mass_max"])
        mask_good_jets = mask_presel & mask_msd & mask_lepton_cleaning

        if jet_tagger != "":
            if "PNetMD" in jet_tagger:
                BB   = "particleNet_XbbVsQCD"
                CC   = "particleNet_XccVsQCD"
            elif "PNet" in jet_tagger:
                BB   = "particleNetWithMass_HbbvsQCD"
                CC   = "particleNetWithMass_HccvsQCD"
            elif "GloParT" in jet_tagger:
                BB = "globalParT3_Xbb"
                CC = "globalParT3_Xcc"
                QCD = "globalParT3_QCD"
            else:
                raise NotImplementedError(f"This tagger is not implemented: {jet_tagger}")
            
            if BB not in jets.fields or CC not in jets.fields:
                raise NotImplementedError(f"{BB} and/or {CC} are not available in the input.")
            if "GloParT" in jet_tagger:
                jets["btagBB"] = jets[BB] / (jets[BB] + jets[QCD])
                jets["btagCC"] = jets[CC] / (jets[CC] + jets[QCD])
            else:
                jets["btagBB"] = jets[BB]
                jets["btagCC"] = jets[CC]

    return jets[mask_good_jets], mask_good_jets

def jet_sel(events, jet_type, params, leptons_collection=""):

    jets = events[jet_type]
    cuts = params.object_preselection["Jet"] #FIXME hardcoded maybe
    # Only jets that are more distant than dr to ALL leptons are tagged as good jets
    # Mask for  jets not passing the preselection
    mask_presel = (
        (jets.pt > cuts["pt"])
        & (np.abs(jets.eta) < cuts["eta"])
        & (jets.jetId >= cuts["jetId"])
    )
    # Lepton cleaning
    if leptons_collection != "":
        dR_jets_lep = jets.metric_table(events[leptons_collection])
        #mask_lepton_cleaning = ak.prod(dR_jets_lep >= cuts["dr_lepton"], axis=2) == 1
        mask_lepton_cleaning = is_lep_cleaned(events, leptons_collection, "Jet", 0.4)

    #if jet_type == "Jet":
        # Selection on PUid. Only available in Run2 UL, thus we need to determine which sample we run over;

    mask_jetpuid = (jets.puId >= 4) | (jets.pt > 50)

    mask_good_jets = mask_presel & mask_lepton_cleaning & mask_jetpuid

    #elif jet_type == "FatJet":
    #    # Apply the msd and preselection cuts
    #    mask_msd = events.FatJet.msoftdrop > cuts["msd"]
    #    mask_good_jets = mask_presel & mask_msd

    return jets[mask_good_jets], mask_good_jets

def fatjet_sel(events, params, leptons_collection=""):

    fatjets = events['FatJet']
    cuts = params.object_preselection['FatJet']

    passes_eta = abs(fatjets.eta) < cuts['eta']
    passes_pt = fatjets.pt > cuts['pt']
    passes_id = fatjets.jetId >= cuts['jetId']
    passes_massLow = fatjets.particleNet_mass >= cuts['mass']['low']
    passes_massHigh = fatjets.particleNet_mass <= cuts['mass']['high']
    if leptons_collection != "":
        dR_jets_lep = fatjets.metric_table(events[leptons_collection])
        #mask_lepton_cleaning = ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1
        mask_lepton_cleaning = is_lep_cleaned(events, leptons_collection, "FatJet", 0.8)

    good_fatjets = passes_eta & passes_pt & passes_id & passes_massLow & passes_massHigh & mask_lepton_cleaning

    return fatjets[good_fatjets]

def fatjet_sel2(events, params, leptons_collection=""):

    fatjets = events['FatJet']
    cuts = params.object_preselection['FatJet']

    passes_eta = abs(fatjets.eta) < cuts['eta']
    passes_pt = fatjets.pt > cuts['pt']
    passes_id = fatjets.jetId >= cuts['jetId']
    passes_massLow = fatjets.particleNet_mass >= cuts['mass']['low']
    passes_massHigh = fatjets.particleNet_mass <= cuts['mass']['high']
    if leptons_collection != "":
        dR_jets_lep = fatjets.metric_table(events[leptons_collection])
        mask_lepton_cleaning = ak.prod(dR_jets_lep > cuts["dr_lepton"], axis=2) == 1
        #mask_lepton_cleaning = is_lep_cleaned(events, leptons_collection, "FatJet", 0.8)

    good_fatjets = passes_eta & passes_pt & passes_id & passes_massLow & passes_massHigh & mask_lepton_cleaning

    return fatjets[good_fatjets]

'''
b-tagging
'''

def bjet_sel(Jet, params):
    return Jet[Jet.btagB > params.object_preselection['btag']['wp']]

'''
q-selection (AK4 below btag WP)
'''

def qjet_sel(Jet, params):
    return Jet[Jet.btagB < params.object_preselection['btag']['wp']]
