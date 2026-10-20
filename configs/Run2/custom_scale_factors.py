import re
import zlib
import numpy as np
import awkward as ak
from pocket_coffea.lib.correction_cache import load_correction_set
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.scale_factors import sf_jet_puId
from omegaconf import OmegaConf


def _load_correction(file, name):
    return load_correction_set(file)[name]


def sf_jet_puId_era(params, jets, year, njets):
    """Pileup jet ID SF with the PU ID bit value of this era's working point.

    pocket-coffea's sf_jet_puId compares the jet puId with the single number
    object_preselection.Jet.puId.value, but the UL bit values differ per era
    (2016: L/M/T = 1/3/7, 2017/2018: 4/6/7). The selection (jet_selection_custom) takes the value
    from jet_scale_factors.jet_puId[year].working_point, so the same value is passed here."""
    puid = OmegaConf.to_container(params.object_preselection["Jet"]["puId"], resolve=True)
    era_cfg = OmegaConf.to_container(params.jet_scale_factors.jet_puId[year], resolve=True)
    puid["value"] = era_cfg["working_point"][puid["wp"]]
    era_params = OmegaConf.create({"object_preselection": {"Jet": {"puId": puid}},
                                   "jet_scale_factors": {"jet_puId": {year: era_cfg}}})
    return sf_jet_puId(era_params, jets, year, njets)

def chunk_rng(events):
    """Random generator seeded by the chunk (file + entry range): reproducible between runs
    and identical for every shape variation of the same chunk."""
    md = events.metadata
    key = f"{md.get('filename', '')}:{md.get('entrystart', 0)}:{md.get('entrystop', 0)}"
    return np.random.default_rng(zlib.crc32(key.encode()))

def btag_jes_systematic(variation):
    """BTV shape-SF systematic for an AK4 JES shape variation of the jet calibrator, e.g.
         'AK4PFchs_JES_TotalUp'                     -> 'up_jes'
         'AK4PFchs_JES_Regrouped_Absolute_2018Down' -> 'down_jesAbsolute_2018'
         'AK4PFchs_JES_FlavorQCDUp'                 -> 'up_jesFlavorQCD'
         'AK4AK8_JES_AbsoluteStatUp'                -> 'up_jesAbsoluteStat'  (AK4+AK8 shifted together)
    None for every other variation (nominal, JER, AK8 jets only, pure-flavour JES, ...), which use
    the central SF."""
    m = re.match(r"^AK4[A-Za-z0-9]*_JES_(.+?)(Up|Down)$", variation)
    if m is None or m.group(1).startswith("FlavorPure"):    # BTV has no pure-flavour JES variations
        return None
    source, direction = m.groups()
    source = source.replace("Regrouped_", "")
    return f"{direction.lower()}_jes{'' if source == 'Total' else source}"


def sf_btag_jes(params, jets, year, njets, systematic):
    """Per-event b-tag shape SF for a JES-varied jet collection: the JES-varied SF ('up_jes*' /
    'down_jes*') for b and light jets, the central SF for c jets (BTV recommendation; the c-jet
    uncertainties are only the cferr ones)."""
    btagSF = params.jet_scale_factors.btagSF[year]
    corr = _load_correction(btagSF.file, btagSF.name)
    discriminator = params.btagging.working_point[year]["btagging_algorithm"]

    flavour = ak.to_numpy(ak.flatten(jets.hadronFlavour))
    abseta = np.abs(ak.to_numpy(ak.flatten(jets.eta)))
    pt = ak.to_numpy(ak.flatten(jets.pt))
    discr = ak.to_numpy(ak.flatten(jets[discriminator]))

    sf = np.asarray(corr.evaluate("central", flavour, abseta, pt, discr), dtype=np.float64)
    notc = flavour != 4
    if np.any(notc):
        try:
            sf[notc] = corr.evaluate(systematic, flavour[notc], abseta[notc], pt[notc], discr[notc])
        except Exception as err:
            raise ValueError(f"b-tag SF systematic '{systematic}' not available in {btagSF.file} "
                             f"({btagSF.name}); check the JES source naming") from err
    return ak.prod(ak.unflatten(sf, njets), axis=1)

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
    

def set_fatjet_mass(fatjets, pnet_mass):
    """Fatjet mass convention of the analysis: every fatjet mass is the ParticleNet regressed mass,
    with JMS/JMR applied in MC (apply_jms_jmr) and as-is in data. Both `mass` (FatJetGood mass
    window, fatjet mass columns, 4-vectors, Xbb SF binning) and `particleNet_mass` (ZH_M,
    outZH_max_ak8pnetMass, mAK8 histogram) carry it, so all of them see the same value."""
    fatjets = ak.with_field(fatjets, pnet_mass, "mass")
    return ak.with_field(fatjets, pnet_mass, "particleNet_mass")


def apply_jms_jmr(events, fatjets, jms_sf, jmr_sf, sig_mc, rng=None):
    """
    Applies JMS and JMR corrections to the ParticleNet regressed mass.
    rng: numpy Generator for the stochastic smearing (see chunk_rng); global numpy RNG if None.
    """
    # 1. Apply Jet Mass Scale (JMS)
    mass_jms = fatjets.particleNet_mass * jms_sf
    
    # 2. Match to GenJetAK8 (Requirement: dR < R/2)
    matched_genjets, matched_fatjets, deltaR = object_matching(
        events.GenJetAK8, fatjets, dr_min=0.4
    )
    
    # Create a mask for valid matches
    has_match = ~ak.is_none(matched_genjets, axis=1)
    gen_mass = ak.fill_none(matched_genjets.mass, 0.0)
    
    # 3. Calculate Deterministic Smearing (Matched)
    mass_deterministic = gen_mass + jmr_sf * (mass_jms - gen_mass)
    
    # 4. Calculate Stochastic Smearing (Unmatched)
    # Generate random numbers from a standard normal distribution N(0,1)
    n_fatjets = int(ak.sum(ak.num(fatjets)))
    gauss = rng.standard_normal(n_fatjets) if rng is not None else np.random.normal(0, 1, n_fatjets)
    random_gauss = ak.unflatten(gauss, ak.num(fatjets))
    
    # Prevent negative arguments in the square root
    smear_factor = np.sqrt(np.maximum(jmr_sf**2 - 1.0, 0.0))
    mass_stochastic = mass_jms * (1.0 + random_gauss * sig_mc * smear_factor)
    
    # 5. Combine using the match mask
    mass_corr = ak.where(has_match, mass_deterministic, mass_stochastic)
    
    # 6. Safety check: prevent negative masses from extreme smearing
    mass_corr = ak.where(mass_corr < 0, 0.0, mass_corr)
    
    return mass_corr


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

# Number of pdfweight_<i> columns written for every MC event (the largest PDF set used, see
# pdf_uncertainty.py); members a sample does not have are NaN
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
    # absolute statistical and systematic parts of the uncertainty (systup/systdown above = total)
    stat = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "stat"
    )
    syst = muon_correctionset[sfName].evaluate(
        np.abs(eta.to_numpy()), pt_.to_numpy(), "syst"
    )
    
    # The unflattened arrays are returned in order to have one row per event.
    return (
        ak.unflatten(sf, counts),
        ak.unflatten(sfup, counts),
        ak.unflatten(sfdown, counts),
        ak.unflatten(stat, counts),
        ak.unflatten(syst, counts),
    )

def sf_mu(params, events, year, key=''):
    '''
    This function computes the per-muon id SF and returns the corresponding per-event SF, obtained by multiplying the per-muon SF in each event.
    Additionally, also the up and down variations of the SF are returned (total uncertainty), and
    a dict with the variations by the statistical and the systematic part alone:
    {'statup', 'statdown', 'systup', 'systdown'} (MUO: syst correlated across years, stat not).
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
    sf, sfup, sfdown, stat, syst = get_mu_sf(params, year, mu_pt_flat, mu_eta_flat, mu_counts, key)

    # The SF arrays corresponding to all the muons are multiplied along the
    # muon axis in order to obtain a per-event scale factor.
    parts = {"statup": sf + stat, "statdown": sf - stat, "systup": sf + syst, "systdown": sf - syst}
    return (ak.prod(sf, axis=1), ak.prod(sfup, axis=1), ak.prod(sfdown, axis=1),
            {name: ak.prod(values, axis=1) for name, values in parts.items()})


