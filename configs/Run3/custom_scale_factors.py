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

# 4FS ttbb (Powheg, Summer24): the nominal of the samples is (muR, muF) = (mu0, 2 mu0) (TMG choice for
# Run 3), but the reweighting grid varies both scales by 0.5 / 1 / 2 around (mu0, mu0) (ttH(bb) report on the
# TTBB theory weights). Following the Run 2 central scale (option 2 of that discussion), the analysis
# nominal of these samples is (mu0, mu0): the event weight gets ttbb_scale_ref = w(mu0, mu0) / w(nominal),
# and the scale variations are taken relative to it. LHEScaleWeight is ordered (muR, muF) in {0.5, 1, 2}^2
# relative to mu0, so (mu0, 2 mu0) = the sample nominal sits at TTBB_SAMPLE_NOMINAL_INDEX (weight 1 by
# construction, checked) and (mu0, mu0) at TTBB_SCALE_REF_INDEX; the standard indices of sf_scaleweights
# then give the 7-point variations around (mu0, mu0) (mu_f_up = the sample nominal).
TTBB_SCALE_FIX_PREFIX = "TTbb_"         # TTbb_Hadronic / _SemiLeptonic / _2L2Nu and their variation samples
TTBB_SCALE_REF_INDEX = 4
TTBB_SAMPLE_NOMINAL_INDEX = 5


def is_ttbb_scale_fix(sample):
    return sample.startswith(TTBB_SCALE_FIX_PREFIX)


def ttbb_scale_ref(events, sample):
    """Per-event factor from the sample nominal to the analysis nominal (mu0, mu0) for the 4FS ttbb samples,
    1 otherwise (numpy). Raises if the scale weights do not have the expected layout."""
    n = len(events)
    if not is_ttbb_scale_fix(sample) or n == 0:
        return np.ones(n)
    if "LHEScaleWeight" not in events.fields:
        raise KeyError(f"{sample}: no LHEScaleWeight, cannot move the ttbb nominal to (mu0, mu0)")
    w = events.LHEScaleWeight
    if not np.all(ak.to_numpy(ak.num(w)) == 9):
        raise ValueError(f"{sample}: expected 9 LHEScaleWeight entries for the ttbb scale fix")
    w = ak.to_numpy(w)
    if not np.allclose(w[:, TTBB_SAMPLE_NOMINAL_INDEX], 1.0, atol=5e-3):
        raise ValueError(f"{sample}: LHEScaleWeight[{TTBB_SAMPLE_NOMINAL_INDEX}] is not the sample nominal "
                         f"(mean {np.mean(w[:, TTBB_SAMPLE_NOMINAL_INDEX]):.4f}): not the scale layout of the "
                         f"remade TTBB samples (check_lhe_weights.py); see custom_scale_factors.TTBB_SCALE_REF_INDEX")
    return w[:, TTBB_SCALE_REF_INDEX]


def scale_ratios(events, sample):
    """sf_scaleweights relative to the analysis nominal: for the 4FS ttbb samples divided by ttbb_scale_ref
    (variations around (mu0, mu0)), otherwise unchanged. None if the sample has no scale weights."""
    scale = sf_scaleweights(events)
    if scale is None or not is_ttbb_scale_fix(sample) or len(events) == 0:
        return scale
    ref = ttbb_scale_ref(events, sample)
    return tuple(ak.to_numpy(x) / ref for x in scale)


# Number of pdfweight_<i> columns written for every MC event (the largest PDF set used, see
# pdf_uncertainty.py); members a sample does not have are NaN. As in Run 2.
N_PDF_COLUMNS = 103

def sf_pdfweights(events):
    """The LHEPdfWeight members (w_member / w_nominal) as a list of numpy arrays, one per member
    (empty list if the sample has no PDF weights). The number of members identifies the PDF set."""
    if "LHEPdfWeight" not in events.fields or len(events) == 0:
        return []
    counts = ak.to_numpy(ak.num(events.LHEPdfWeight))
    if counts.min() != counts.max():
        raise ValueError(f"LHEPdfWeight has a varying number of members ({counts.min()}-{counts.max()})")
    pdf = ak.to_numpy(events.LHEPdfWeight)
    return [pdf[:, i] for i in range(counts[0])]

def pdfweight_columns(events, n_columns=N_PDF_COLUMNS):
    """{'pdfweight_<i>': array} for i < n_columns, NaN where the sample has fewer members."""
    members = sf_pdfweights(events)
    if len(members) > n_columns:
        raise ValueError(f"{len(members)} LHEPdfWeight members, only {n_columns} columns")
    missing = np.full(len(events), np.nan)
    return {f"pdfweight_{i}": members[i] if i < len(members) else missing for i in range(n_columns)}

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
    (data and MC; the JMS/JMR variations are applied on top in the processor), used by the FatJetGood
    mass window, ZH_M, the fatjet mass columns and 4-vectors. The JEC-corrected jet mass is kept as
    FatJet.mass_jec."""
    if "mass_jec" not in fatjets.fields:
        fatjets = ak.with_field(fatjets, fatjets.mass, "mass_jec")
    return ak.with_field(fatjets, gloparT_regressed_mass(fatjets), "mass")


# ---------------------------------------------------------------- b tagging (fixed WPs)
# Uncertainty sources of the multi-WP b-tag SFs: column suffix -> (jets varied, BTV systematic).
# b and c jets (UParTAK4_comb, c = b SF with doubled syst) vary together; light jets (UParTAK4_light)
# separately. 'correlated' (syst) is correlated between years, 'uncorrelated' (stat) is not.
BTAG_SF_SOURCES = {
    "bc_correlated": ("bc", "correlated"),
    "bc_uncorrelated": ("bc", "uncorrelated"),
    "light_correlated": ("light", "correlated"),
    "light_uncorrelated": ("light", "uncorrelated"),
}


def btag_sf_columns():
    """Output columns of sf_btag_multiwp: btag_sf and btag_sf_<source>_up/down."""
    return ["btag_sf"] + [f"btag_sf_{src}_{d}" for src in BTAG_SF_SOURCES for d in ("up", "down")]


def sf_btag_multiwp(params, jets, year, wps, eff_file, eff_group):
    """Fixed-WP b-tag SFs, BTV method 1a extended to several working points:

        w = prod_{tagged at the tightest WP} SF_n
          * prod_{tagged at WP k, not at WP k+1} (SF_k eps_k - SF_k+1 eps_k+1) / (eps_k - eps_k+1)
          * prod_{not tagged at the loosest WP}   (1 - SF_1 eps_1) / (1 - eps_1)

    over all `jets` (the jets whose b-tag information the analysis uses), with eps the MC
    efficiencies (eff_file: correction UParTAK4_eff_<eff_group>, from derive_btag_eff.py) and SF
    from BTV (UParTAK4_comb for b/c, UParTAK4_light for light jets). A factor with a negative
    numerator (SF_k+1 eps_k+1 > SF_k eps_k, unphysical) is set to 1, as BTV recommends; the same is
    done for the untagged factor (SF_1 eps_1 > 1) and for a factor with a vanishing denominator.

    wps: working points, loosest first (analysis_years.BTAG_SF_WPS). A jet's tag level is the number of
    them it passes (btagB > cut), so it agrees with JetGood.btag_<WP> and bJetGood.
    Returns {column: per-event weight} for btag_sf_columns().
    """
    cfg = params.jet_scale_factors.btagSF[year]
    sf_file, tagger = cfg["file"], cfg["name"]
    cuts = params.btagging.working_point[year]["btagging_WP"]
    if any(cuts[a] >= cuts[b] for a, b in zip(wps, wps[1:])):
        raise ValueError(f"b-tag working points {wps} not ordered loosest first: {dict(cuts)}")
    sf_bc = _load_correction(sf_file, f"{tagger}_comb")
    sf_light = _load_correction(sf_file, f"{tagger}_light")
    eff = _load_correction(eff_file, f"{tagger}_eff_{eff_group}")

    counts = ak.num(jets)
    flat = ak.flatten(jets)
    pt = ak.to_numpy(flat.pt).astype(np.float64)
    abseta = np.abs(ak.to_numpy(flat.eta)).astype(np.float64)
    flavor = ak.to_numpy(flat.hadronFlavour).astype(int)
    btag = ak.to_numpy(flat.btagB)
    level = sum((btag > cuts[wp]).astype(int) for wp in wps)       # 0 = untagged, len(wps) = tightest
    heavy = flavor > 0

    eps = {wp: np.clip(eff.evaluate(wp, flavor, abseta, pt), 0.0, 1.0) for wp in wps}

    def scale_factors(syst_bc, syst_light):
        out = {}
        for wp in wps:
            sf = np.ones(len(pt))
            if heavy.any():
                sf[heavy] = sf_bc.evaluate(syst_bc, wp, flavor[heavy], abseta[heavy], pt[heavy])
            if (~heavy).any():
                sf[~heavy] = sf_light.evaluate(syst_light, wp, flavor[~heavy], abseta[~heavy], pt[~heavy])
            out[wp] = sf
        return out

    def event_weight(sf):
        num, den = 1.0 - sf[wps[0]] * eps[wps[0]], 1.0 - eps[wps[0]]
        factor = np.where((num > 0) & (den > 0), num / np.where(den > 0, den, 1.0), 1.0)
        factor = np.where(level == 0, factor, 1.0)
        for k in range(len(wps) - 1):
            lo, hi = wps[k], wps[k + 1]
            num = sf[lo] * eps[lo] - sf[hi] * eps[hi]
            den = eps[lo] - eps[hi]
            f = np.where((num > 0) & (den > 0), num / np.where(den > 0, den, 1.0), 1.0)
            factor = np.where(level == k + 1, f, factor)
        factor = np.where(level == len(wps), sf[wps[-1]], factor)
        return ak.to_numpy(ak.prod(ak.unflatten(factor, counts), axis=1))

    out = {"btag_sf": event_weight(scale_factors("central", "central"))}
    for src, (jets_varied, kind) in BTAG_SF_SOURCES.items():
        for d in ("up", "down"):
            syst = f"{d}_{kind}"
            sf = scale_factors(syst, "central") if jets_varied == "bc" else scale_factors("central", syst)
            out[f"btag_sf_{src}_{d}"] = event_weight(sf)
    return out
