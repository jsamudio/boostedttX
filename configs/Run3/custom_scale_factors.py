"""Scale factors and weights of the Run 3 analysis."""
import awkward as ak
import numpy as np
from pocket_coffea.lib.correction_cache import load_correction_set


def _load_correction(file, name):
    return load_correction_set(file)[name]

def sf_lep_trig(events, file, collection, correction, variations=("nominal", "up", "down")):
    evaluator = _load_correction(file, correction)
    lep = events[collection]
    counts = ak.num(lep.pt)
    pt = ak.to_numpy(ak.flatten(lep.pt)).astype(np.float64)
    eta = ak.to_numpy(ak.flatten(lep.eta)).astype(np.float64)
    sfs = [ak.unflatten(evaluator.evaluate(pt, eta, v), counts) for v in variations]
    return tuple(ak.prod(sf, axis=1) for sf in sfs)

def sf_ele_trig(events, file, **kw):
    return sf_lep_trig(events, file, "ElectronGood", "EleTrigSF", **kw)

def sf_mu_trig(events, file, **kw):
    return sf_lep_trig(events, file, "MuonGood", "MuTrigSF", **kw)

def sf_scaleweights(events):
    # Ensure this is MC, not Data
    if not hasattr(events, "LHEScaleWeight") or events.LHEScaleWeight is None:
        return None # Handle data appropriately
    
    n_weights = ak.num(events.LHEScaleWeight)[0]
    
    if n_weights == 9:
        mu_r_up    = events.LHEScaleWeight[:, 7]
        mu_r_down  = events.LHEScaleWeight[:, 1]
        mu_f_up    = events.LHEScaleWeight[:, 5]
        mu_f_down  = events.LHEScaleWeight[:, 3]
        mu_rf_up   = events.LHEScaleWeight[:, 8]
        mu_rf_down = events.LHEScaleWeight[:, 0]
    elif n_weights == 8:
        # Nominal is missing, indices > 3 shift down by 1
        mu_r_up    = events.LHEScaleWeight[:, 6]
        mu_r_down  = events.LHEScaleWeight[:, 1]
        mu_f_up    = events.LHEScaleWeight[:, 4]
        mu_f_down  = events.LHEScaleWeight[:, 3]
        mu_rf_up   = events.LHEScaleWeight[:, 7]
        mu_rf_down = events.LHEScaleWeight[:, 0]
    else:
        raise ValueError(f"Unexpected number of LHEScaleWeights: {n_weights}")

    return mu_r_up, mu_r_down, mu_f_up, mu_f_down, mu_rf_up, mu_rf_down

def sf_partonshower_isr(events):
    '''Up and down variations for the ISR parton shower weights.
    In order to properly store the weights, a dummy weight of 1 is stored
    as central value for the ISR correction.
    Conventions for the PS weights are:
    [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
    '''
    isr_down = events.PSWeight[:,2]
    isr_up = events.PSWeight[:,0]
    nom = ak.ones_like(isr_up)

    return nom, isr_up, isr_down

def sf_partonshower_fsr(events):
    '''Up and down variations for the FSR parton shower weights.
    In order to properly store the weights, a dummy weight of 1 is stored
    as central value for the FSR correction.
    Conventions for the PS weights are:
    [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
    '''
    fsr_down = events.PSWeight[:,3]
    fsr_up = events.PSWeight[:,1]
    nom = ak.ones_like(fsr_up)

    return nom, fsr_up, fsr_down

def get_mu_sf(params, year, pt, eta, counts, key=''):
    '''
    This function computes the per-muon id or iso SF.
    '''
    muonSF = params["lepton_scale_factors"]["muon_sf"]

    muon_correctionset = load_correction_set(muonSF.JSONfiles[year]['file'])

    if key not in ["id","iso","trigger", "reco"]:
        raise Exception(f"Muon SF key {key} not recognized")
    
    sfName = muonSF.sf_name[year][key]

    pt_ = pt

    if key == 'reco':
        pt_ = ak.where(pt_ < 40, 40, pt_) # for RECO SF, binning starts at 40 GeV, and recommendation is to use the 40 bin for anything between 10 and 40 GeV
    
    sf = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "nominal"
    )
    sfup = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "systup"
    )
    sfdown = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "systdown"
    )
    stat = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "systdown"
    )
    
    # The unflattened arrays are returned in order to have one row per event.
    return (
        ak.unflatten(sf, counts),
        ak.unflatten(sfup, counts),
        ak.unflatten(sfdown, counts),
        ak.unflatten(stat, counts)
    )

def sf_mu(params, events, year, key=''):
    '''
    This function computes the per-muon id SF and returns the corresponding per-event SF, obtained by multiplying the per-muon SF in each event.
    Additionally, also the up and down variations of the SF are returned.
    '''
    coll = params.lepton_scale_factors.muon_sf.collection
    mu_pt = events[coll].pt
    mu_eta = events[coll].eta

    # Since `correctionlib` does not support jagged arrays as an input, the pt and eta arrays are flattened.
    mu_pt_flat, mu_eta_flat, mu_counts = (
        ak.flatten(mu_pt),
        ak.flatten(mu_eta),
        ak.num(mu_pt),
    )
    sf, sfup, sfdown, stat = get_mu_sf(params, year, mu_pt_flat, mu_eta_flat, mu_counts, key)

    # The SF arrays corresponding to all the muons are multiplied along the
    # muon axis in order to obtain a per-event scale factor.
    return ak.prod(sf, axis=1), ak.prod(sfup, axis=1), ak.prod(sfdown, axis=1), ak.prod(stat, axis=1)


def gloparT_regressed_mass(fatjets):
    """GloParT 2-prong regressed mass: correction factor times the raw (uncorrected) jet mass.
    mass * (1 - rawFactor) is the raw mass also after pocket-coffea's JEC re-application.
    Starts from the JEC mass saved by set_fatjet_mass if present, so it never compounds."""
    jec_mass = fatjets.mass_jec if "mass_jec" in fatjets.fields else fatjets.mass
    return fatjets.globalParT3_massCorrX2p * jec_mass * (1.0 - fatjets.rawFactor)


def set_fatjet_mass(fatjets):
    """Fatjet mass convention of the analysis: FatJet.mass is the GloParT 2-prong regressed mass
    (data and MC; no JMS/JMR for it yet), used by the FatJetGood mass window, ZH_M, the fatjet
    mass columns and 4-vectors. The JEC-corrected jet mass is kept as FatJet.mass_jec."""
    if "mass_jec" not in fatjets.fields:
        fatjets = ak.with_field(fatjets, fatjets.mass, "mass_jec")
    return ak.with_field(fatjets, gloparT_regressed_mass(fatjets), "mass")
