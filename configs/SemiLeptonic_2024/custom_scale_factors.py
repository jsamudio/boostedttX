import numpy as np
import awkward as ak
import correctionlib
from coffea.lookup_tools import txt_converters, rochester_lookup

def sf_ele_trig(events):
    # Load the correction set once at the start of your processor
    cset = correctionlib.CorrectionSet.from_file("/cms/data/jsamudio/boosted/boostedttX/configs/SemiLeptonic_2024/ele_trig_sf_2024.json.gz")
    evaluator = cset["EleTrigSF"]

    # In your event loop (assuming you have flat awkward arrays of lepton pT and Eta)
    ele_pt = events.ElectronGood.pt
    ele_eta = events.ElectronGood.eta

    ele_pt_flat, ele_eta_flat, ele_counts = (
        ak.flatten(ele_pt).to_numpy(),
        ak.flatten(ele_eta).to_numpy(),
        ak.num(ele_pt),
    )

    scale_factors = [
        ak.unflatten(
            evaluator.evaluate(ele_eta_flat, ele_pt_flat, variation),
            ele_counts,
        )
        for variation in ("nominal", "up", "down")
    ]
    

    # return a per-event scale factor by multiplying all electron scale factors
    return tuple(ak.prod(sf, axis=1) for sf in scale_factors)


    

def sf_bbtag(params, fatjets, year, njets, variations=["central"]):
    '''
    Application of custom bbtag SFs. These should be applied to fatjets passing the WP threshold.
    '''
    bbtagSF = params.jet_scale_factors.bbtagSF[year]
    working_point = bbtagSF.wp
    #bbtag_discriminator = params.bbtagging.working_point[year]["bbtagging_algorithm"]
    cset = correctionlib.CorrectionSet.from_file(bbtagSF.file)
    corr = cset[bbtagSF.name]

    #flavour = ak.to_numpy(ak.flatten(jets.hadronFlavour))
    #abseta = np.abs(ak.to_numpy(ak.flatten(jets.eta)))
    pt = ak.to_numpy(ak.flatten(fatjets.pt))
    score = ak.to_numpy(ak.flatten(fatjets.xbbVsQCD))
    #discr = ak.to_numpy(ak.flatten(jets[bbtag_discriminator]))

    central_SF_byjet = corr.evaluate("central", working_point, pt)

    def _get_sf_variation_with_mask(variation, mask):
        index = (np.indices(score.shape)).flatten()[mask]
        # Copying the central SF
        sf = np.copy(central_SF_byjet)
        w = corr.evaluate(variation, working_point, pt[mask])
        sf[index] = w
        sf_out = ak.prod(ak.unflatten(sf, njets), axis=1)
        return sf_out

    output = {}
    for variation in variations:
        # FIXME hardcoded score cutoff
        # but this should make it so the SF only applies to real signal
        score_mask = score > 0.9105
        output[variation] = [_get_sf_variation_with_mask(f"{variation}", score_mask)]
    print(output)

    return output

def sf_scaleweights(events):
    '''Up and down variations for the ISR parton shower weights.
    In order to properly store the weights, a dummy weight of 1 is stored
    as central value for the ISR correction.
    Conventions for the PS weights are:
    [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
    '''
    mu_r_up    = events.LHEScaleWeight[:,7]
    mu_r_down  = events.LHEScaleWeight[:,1]
    mu_f_up    = events.LHEScaleWeight[:,5]
    mu_f_down  = events.LHEScaleWeight[:,3]
    mu_rf_up   = events.LHEScaleWeight[:,8]
    mu_rf_down = events.LHEScaleWeight[:,0]

    return mu_r_up, mu_r_down, mu_f_up, mu_f_down, mu_rf_up, mu_rf_down,

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

    muon_correctionset = correctionlib.CorrectionSet.from_file(
        muonSF.JSONfiles[year]['file']
    )

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

def ApplyRochesterCorrections(year, mu, is_data):
    if year.startswith('201'): #Run2 scenario
        rocco_tag = None
        if year == '2016':
            rocco_tag = "2016bUL"
        elif year == '2016APV':
            rocco_tag = "2016aUL"
        elif year == '2017':
            rocco_tag = "2017UL"
        elif year == '2018':
            rocco_tag = "2018UL"
        rochester_data = txt_converters.convert_rochester_file(f"/cms/data/jsamudio/boosted/boostedttX/RoccoR/RoccoR{rocco_tag}.txt", loaduncs=True) # FIXME get the right file path when imported
        rochester = rochester_lookup.rochester_lookup(rochester_data)
        if not is_data:
            hasgen = ~np.isnan(ak.fill_none(mu.matched_gen.pt, np.nan))
            mc_rand = np.random.rand(*ak.to_numpy(ak.flatten(mu.pt)).shape)
            mc_rand = ak.unflatten(mc_rand, ak.num(mu.pt, axis=1))
            corrections = np.array(ak.flatten(ak.ones_like(mu.pt)))
            corrections_unc = np.array(ak.flatten(ak.ones_like(mu.pt)))
            mc_kspread = rochester.kSpreadMC(
                mu.charge[hasgen],mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen]
            )
            mc_kspread_stat = ak.std(ak.concatenate([rochester.kSpreadMC(
                mu.charge[hasgen],
                mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen],
                s=1,
                m=i
            ) for i in range(100)],axis=-1), axis=-1)
            mc_kspread_zpt = abs(mc_kspread - rochester.kSpreadMC(
                mu.charge[hasgen],
                mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen],
                s=2
            ))
            mc_kspread_ewk = abs(mc_kspread - rochester.kSpreadMC(
                mu.charge[hasgen],
                mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen],
                s=3
            ))
            mc_kspread_deltam = abs(mc_kspread - rochester.kSpreadMC(
                mu.charge[hasgen],
                mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen],
                s=4
            ))
            mc_kspread_ewk2 = abs(mc_kspread - rochester.kSpreadMC(
                mu.charge[hasgen],
                mu.pt[hasgen],
                mu.eta[hasgen],
                mu.phi[hasgen],
                mu.matched_gen.pt[hasgen],
                s=5
            ))
            mc_kspread_unc = np.sqrt(mc_kspread_stat**2 + mc_kspread_ewk**2 + mc_kspread_zpt**2 + mc_kspread_ewk2**2 + mc_kspread_deltam**2)
            #print("Spread unc:", mc_kspread_stat**2)
            mc_ksmear = rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen]
            )
            mc_ksmear_stat = ak.std(ak.concatenate([rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
                s=1,
                m=i
            ) for i in range(100)],axis=-1), axis=-1)
            mc_ksmear_zpt = abs(mc_ksmear - rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
                s=2
            ))
            mc_ksmear_ewk = abs(mc_ksmear - rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
                s=3
            ))
            mc_ksmear_deltam = abs(mc_ksmear - rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
                s=4
            ))
            mc_ksmear_ewk2 = abs(mc_ksmear - rochester.kSmearMC(
                mu.charge[~hasgen],
                mu.pt[~hasgen],
                mu.eta[~hasgen],
                mu.phi[~hasgen],
                mu.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
                s=5
            ))
            mc_ksmear_unc = np.sqrt(mc_ksmear_stat**2 + mc_ksmear_ewk**2 + mc_ksmear_zpt**2 + mc_ksmear_ewk2**2 + mc_ksmear_deltam**2)
            hasgen_flat = np.array(ak.flatten(hasgen))
            corrections_unc[hasgen_flat] = np.array(ak.flatten(mc_kspread_unc))
            corrections_unc[~hasgen_flat] = np.array(ak.flatten(mc_ksmear_unc))
            corrections[hasgen_flat] = np.array(ak.flatten(mc_kspread))
            corrections[~hasgen_flat] = np.array(ak.flatten(mc_ksmear))
            corrections = ak.unflatten(corrections, ak.num(mu.pt, axis=1))
            corrections_unc = ak.unflatten(corrections_unc, ak.num(mu.pt, axis=1))
        else:
            corrections = rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)
            data_scale_stat = ak.std(ak.concatenate([rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi, s=1, m=i) for i in range(100)], axis=-1), axis=-1)
            data_scale_zpt = abs(corrections - rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi, s=2))
            data_scale_ewk = abs(corrections - rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi, s=3))
            data_scale_deltam = abs(corrections - rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi, s=4))
            data_scale_ewk2 = abs(corrections - rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi, s=5))
            data_unc = np.sqrt(data_scale_stat**2 + data_scale_ewk**2 + data_scale_zpt**2 + data_scale_deltam**2 + data_scale_ewk2**2)
            
    else:
        corrections = ak.ones_like(mu.pt)
    if not is_data:
        return (mu.pt * corrections), corrections_unc # FIXME double check the output but when it goes to 
    else:
        return (mu.pt * corrections), ak.ones_like(mu.pt)

JECjsonFiles = {
    '2016_PreVFP': {
        'AK4': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016preVFP_UL/jet_jerc.json.gz',
        'AK8': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016preVFP_UL/fatJet_jerc.json.gz',
    },
    '2016_PostVFP': {
        'AK4': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016postVFP_UL/jet_jerc.json.gz',
        'AK8': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016postVFP_UL/fatJet_jerc.json.gz',
    },
    '2017': {
        'AK4': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2017_UL/jet_jerc.json.gz',
        'AK8': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2017_UL/fatJet_jerc.json.gz',
    },
    '2018': {
        'AK4': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2018_UL/jet_jerc.json.gz',
        'AK8': '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2018_UL/fatJet_jerc.json.gz',
    },
}

jec_variations = ['AbsoluteMPFBias', 'AbsoluteScale', 'AbsoluteStat',
                  'FlavorQCD', 'Fragmentation', 'PileUpDataMC',
                  'PileUpPtBB', 'PileUpPtEC1', 'PileUpPtEC2',
                  'PileUpPtHF', 'PileUpPtRef', 'RelativeFSR',
                  'RelativeJEREC1', 'RelativeJEREC2', 'RelativeJERHF',
                  'RelativePtBB', 'RelativePtEC1', 'RelativePtEC2',
                  'RelativePtHF', 'RelativeBal', 'RelativeSample',
                  'RelativeStatEC', 'RelativeStatFSR', 'RelativeStatHF',
                  'SinglePionECAL', 'SinglePionHCAL', 'TimePtEta']

def jet_correction_correctionlib(
    events,
    Jet,
    typeJet,
    year,
    JECversion,
    JERversion=None,
    verbose=False,
    pt_variation="",
    jer_variation="",
):
    """
    Implements Jet Energy Corrections (L1L2L3) and Jet Energy Smearing (Hybrid Method).
    Uses correctionlib and common-POG JSON files.
    """
    
    # --- 1. Setup JEC File ---
    # Determine the correct key for JECjsonFiles based on typeJet (AK4 vs AK8)
    jet_type_key = [t for t in ["AK4", "AK8"] if typeJet.startswith(t)][0]
    jsonfile = JECjsonFiles[year][jet_type_key]
    print(jsonfile)
    
    cset = correctionlib.CorrectionSet.from_file(jsonfile)
    corr = cset.compound[f"{JECversion}_L1L2L3Res_{typeJet}"]

    # --- 2. Prepare Jet Data ---
    jets = events[Jet]
    
    # Calculate raw values to undo previous corrections if necessary
    jets["pt_raw"] = (1 - jets["rawFactor"]) * jets["pt"]
    jets["mass_raw"] = (1 - jets["rawFactor"]) * jets["mass"]
    
    # Broadcast fixedGridRhoFastjetAll to match jet structure
    jets["rho"] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, jets.pt)[0]

    # Flatten for correctionlib (it requires flat arrays)
    j_flat, n_jets = ak.flatten(jets), ak.num(jets)

    # --- 3. Apply Nominal JEC ---
    flat_corr_factor = corr.evaluate(
        np.array(j_flat["area"]),
        np.array(j_flat["eta"]),
        np.array(j_flat["pt_raw"]),
        np.array(j_flat["rho"]),
    )
    
    corr_factor = ak.unflatten(flat_corr_factor, n_jets)

    jets_corrected = ak.copy(jets)
    jets_corrected["pt"] = jets["pt_raw"] * corr_factor
    jets_corrected["mass"] = jets["mass_raw"] * corr_factor
    jets_corrected["rho"] = jets["rho"]

    # --- 4. Apply JEC Uncertainty (if requested) ---
    if pt_variation:
        # Clean the variation string (remove _up/_down suffixes)
        base_variation = pt_variation.replace("_up", "").replace("_down", "")
        corr_unc_provider = cset[f"{JECversion}_{base_variation}_{typeJet}"]
        
        # Evaluate uncertainty using corrected pT
        flat_unc = corr_unc_provider.evaluate(
            np.array(j_flat["eta"]),
            np.array(j_flat["pt_raw"]) * flat_corr_factor,
        )
        
        unc_factor = ak.unflatten(flat_unc, n_jets)
        
        if "up" in pt_variation:
            shift = 1.0 + unc_factor
        else:
            shift = 1.0 - unc_factor
            
        jets_corrected["pt"] = jets_corrected["pt"] * shift
        jets_corrected["mass"] = jets_corrected["mass"] * shift

    # Helper for verbose printing
    seed = events.event[0]
    if verbose:
        print(f"\n{seed} JEC: Untransformed pt ratios:\n{jets.pt[0] / jets.pt_raw[0]}")
        print(f"{seed} JEC: Corrected pt ratios:\n{jets_corrected.pt[0] / jets_corrected.pt_raw[0]}\n")

    # --- 5. Apply JER Smearing (if requested) ---
    if not JERversion:
        return jets_corrected, {}

    # Setup JER tools
    sf_provider = cset[f"{JERversion}_ScaleFactor_{typeJet}"]
    res_provider = cset[f"{JERversion}_PtResolution_{typeJet}"]

    # Flatten corrected jets for evaluation
    j_corr_flat = ak.flatten(jets_corrected)
    
    # 5a. Calculate Scale Factors (SF)
    sf_flat_nom = sf_provider.evaluate(j_corr_flat["eta"].to_numpy(), "nom")
    sf_flat_final = sf_flat_nom  # Default to nominal

    # Handle Split JER variations
    if jer_variation and not pt_variation:
        var_tag = "up" if "up" in jer_variation else "down"
        
        eta_abs = np.abs(j_corr_flat["eta"].to_numpy())
        pt_val = np.abs(j_corr_flat["pt"].to_numpy())
        
        # Define the Split JER regions
        mask = None
        if "1" in jer_variation:
            mask = eta_abs < 1.93
        elif "2" in jer_variation:
            mask = (eta_abs >= 1.93) & (eta_abs < 2.5)
        elif "3" in jer_variation:
            mask = (eta_abs >= 2.5) & (eta_abs < 3) & (pt_val < 50)
        elif "4" in jer_variation:
            mask = (eta_abs >= 2.5) & (eta_abs < 3) & (pt_val >= 50)
        elif "5" in jer_variation:
            mask = (eta_abs >= 3) & (eta_abs < 5) & (pt_val < 50)
        elif "6" in jer_variation:
            mask = (eta_abs >= 3) & (eta_abs < 5) & (pt_val >= 50)

        if mask is not None:
            sf_flat_var = sf_provider.evaluate(j_corr_flat["eta"].to_numpy(), var_tag)
            sf_flat_final = np.where(mask, sf_flat_var, sf_flat_nom)

    scale_factor = ak.unflatten(sf_flat_final, n_jets)

    # 5b. Calculate Resolution
    res_flat = res_provider.evaluate(
        j_corr_flat["eta"].to_numpy(), 
        j_corr_flat["pt"].to_numpy(), 
        j_corr_flat["rho"].to_numpy()
    )
    pt_resolution = ak.unflatten(res_flat, n_jets)

    # 5c. Gen Matching
    # Constants for matching
    dr_min = 0.2 if "AK4" in typeJet else 0.4
    gen_jet_coll = "GenJet" if "AK4" in typeJet else "GenJetAK8"
    gen_jet_idx = "genJetIdx" if "AK4" in typeJet else "genJetAK8Idx"

    # Fetch GenJets
    genjets = events[gen_jet_coll]
    
    # Mask invalid indices (-1 or out of bounds)
    valid_gen_idx = (jets_corrected[gen_jet_idx] >= 0) & (jets_corrected[gen_jet_idx] < ak.num(genjets))
    matched_gen_indices = ak.mask(jets_corrected[gen_jet_idx], valid_gen_idx)
    
    # Fetch the actual matched objects
    matched_genjets = genjets[matched_gen_indices]
    matched_jets = ak.mask(jets_corrected, ~ak.is_none(matched_genjets, axis=1))

    # Apply dPt Check (3 sigma)
    dPt = np.abs(matched_jets.pt - matched_genjets.pt)
    pt_min_req = 3 * pt_resolution * jets_corrected["pt"]
    
    # Final masks for valid matches
    #is_matched = ((~ak.is_none(matched_jets.pt, axis=1)) & (dPt < pt_min_req)).fill_none(False)

    # 1. Create the mask (this will contain True, False, and None)
    is_matched_mask = (~ak.is_none(matched_jets.pt, axis=1)) & (dPt < pt_min_req)
    
    # 2. Use ak.fill_none function to replace None with False
    is_matched = ak.fill_none(is_matched_mask, False)
    
    # 5d. Calculate Smear Factors
    # 1. Deterministic Smearing (Scaling Method)
    det_smear = 1 + (scale_factor - 1) * (matched_jets["pt"] - matched_genjets["pt"]) / matched_jets["pt"]
    
    # 2. Stochastic Smearing
    # Create deterministic seed based on chunk info
    seed_dict = {}
    if "filename" in events.metadata:
        fname = events.metadata["filename"]
        estart = events.metadata.get("entrystart", 0)
        estop = events.metadata.get("entrystop", 0)
        seed_dict[f"chunk_{fname}_{estart}-{estop}"] = seed

    np.random.seed(seed)
    rand_gauss = np.random.normal(np.zeros_like(res_flat), res_flat)
    jer_smear_noise = ak.unflatten(rand_gauss, n_jets)
    
    sqrt_arg = np.maximum(scale_factor**2 - 1, 0)
    stoch_smear = 1 + jer_smear_noise * np.sqrt(sqrt_arg)

    # Combine methods
    final_smear_factor = ak.where(is_matched, det_smear, stoch_smear)
    
    # Apply to jets
    jets_smeared = ak.copy(jets_corrected)
    jets_smeared["pt"] = jets_corrected["pt"] * final_smear_factor
    jets_smeared["mass"] = jets_corrected["mass"] * final_smear_factor

    if verbose:
        print(f"{seed} JER: Smear Factors:\n{final_smear_factor}")
        print(f"{seed} JER: Final Smeared pt ratios:\n{jets_smeared.pt / jets_corrected.pt}\n")

    return jets_smeared, seed_dict

def recompute_type1_met_correctionlib(
    events,
    jets_corrected, # The output from your JEC function
    year,           # The correctionlib set
    JECversion,
    jet_type="AK4PFchs"
):
    """
    Fully recomputes Type-1 MET matching the CMS recipe:
    MET_Final = RawMET - Sum(Jet_Corrected - Jet_L1)
    """

    # --- 1. Prepare Inputs ---
    # We need RawMET, not the standard MET
    met_raw = events.RawMET
    
    # In a full analysis, you must merge 'Jet' and 'CorrT1METJet' here.
    # For simplicity, if you only have 'Jet' corrected, we proceed with that,
    # but acknowledge this is an approximation compared to the full snippet.
    jets = jets_corrected 
    
    # --- 2. Calculate L1 Corrections ---
    # We need to isolate the L1 factor to subtract it. 
    # Type-1 MET = RawMET + (RawJet * L1 - RawJet * L1*L2*L3...) 
    # Wait, the formula is: MET_Type1 = MET_Raw - (Jet_Full - Jet_L1)
    
    # We need to evaluate just the L1 correction for these jets.
    # Assuming 'jets' has 'pt_raw', 'area', 'eta', 'rho' from your previous function.
    jet_type_key = [t for t in ["AK4", "AK8"] if jet_type.startswith(t)][0]
    jsonfile = JECjsonFiles[year][jet_type_key]
    
    cset = correctionlib.CorrectionSet.from_file(jsonfile)
    print(jsonfile)
    l1_corr_provider = cset[f"{JECversion}_L1FastJet_{jet_type}"]
    
    j_flat = ak.flatten(jets)
    flat_l1_factor = l1_corr_provider.evaluate(
        np.array(j_flat["area"]),
        np.array(j_flat["eta"]),
        np.array(j_flat["pt_raw"]),
        np.array(j_flat["rho"]),
    )
    l1_factor = ak.unflatten(flat_l1_factor, ak.num(jets))
    
    # Calculate vector magnitudes
    pt_l1 = jets["pt_raw"] * l1_factor
    pt_final = jets["pt"] # This comes from your main JEC function (L1L2L3+Res)

    # --- 3. Apply MET Cleaning Cuts (Match your snippet) ---
    # pT > 15, |eta| < 5.2, EM Fraction < 0.9
    # Note: NanoAOD branch names for EM frac might differ slightly (neEmEF, chEmEF)
    pass_clean = (
        (pt_final > 15.0) & 
        (abs(jets.eta) < 5.2) & 
        ((jets.neEmEF + jets.chEmEF) < 0.9)
    )

    # --- 4. Calculate Vector Shift ---
    # Shift = Final - L1
    # We only sum over jets that pass the cleaning cuts
    pt_diff = ak.where(pass_clean, pt_final - pt_l1, 0.0)
    
    shift_px = ak.sum(pt_diff * np.cos(jets.phi), axis=1)
    shift_py = ak.sum(pt_diff * np.sin(jets.phi), axis=1)

    # --- 5. Apply to Raw MET ---
    # MET_new = MET_raw - Sum(Vector_Diff)
    new_met_px = met_raw.phi * np.cos(met_raw.phi) # Wait, RawMET usually has pt/phi
    new_met_px = met_raw.pt * np.cos(met_raw.phi) - shift_px
    new_met_py = met_raw.pt * np.sin(met_raw.phi) - shift_py

    new_met_pt = np.hypot(new_met_px, new_met_py)
    new_met_phi = np.arctan2(new_met_py, new_met_px)

    return ak.zip(
        {"pt": new_met_pt, "phi": new_met_phi},
        with_name="PtEtaPhiMLorentzVector"
    )

import numpy as np
import awkward as ak
import correctionlib

import numpy as np
import awkward as ak
import correctionlib

def sf_btag_wp(params, jets, year, njets, mc_efficiencies, wp_threshold, variations=["central"], working_point='M'):
    btagSF = params.jet_scale_factors.btagSF[year]
    btag_discriminator = params.btagging.working_point[year]["btagging_algorithm"]
    cset = correctionlib.CorrectionSet.from_file(btagSF.file)
    
    # 1. LOAD BOTH CORRECTION OBJECTS
    corr_bc = cset["UParTAK4_comb"]   # Heavy flavor (b, c)
    corr_light = cset["UParTAK4_light"]  # Light flavor (udsg)

    # Flatten inputs
    flavour = ak.to_numpy(ak.flatten(jets.hadronFlavour))
    abseta = np.abs(ak.to_numpy(ak.flatten(jets.eta)))
    pt = ak.to_numpy(ak.flatten(jets.pt))
    discr = ak.to_numpy(ak.flatten(jets[btag_discriminator]))
    
    eff = ak.to_numpy(ak.flatten(mc_efficiencies))

    # Evaluate WP Pass/Fail Mask
    is_tagged = discr >= wp_threshold

    def _get_method1a_weight(variation_str, mask):
        # Force any weird flavors to 0
        clean_flavour = np.where((flavour != 5) & (flavour != 4), 0, flavour)
        
        # Reconstruct the full SF array for this variation, default to 1.0
        full_sfs = np.ones_like(flavour, dtype=float)
        
        if np.any(mask):
            # Extract the subset of jets passing the mask
            f_sub = clean_flavour[mask]
            a_sub = abseta[mask]
            p_sub = pt[mask]
            
            # Prepare an array to hold the raw SFs for this subset
            raw_sfs = np.ones_like(f_sub, dtype=float)
            
            # Create sub-masks for heavy vs light flavor
            bc_idx = (f_sub == 5) | (f_sub == 4)
            l_idx = (f_sub == 0)
            
            # EVALUATE HEAVY AND LIGHT SEPARATELY
            if np.any(bc_idx):
                raw_sfs[bc_idx] = corr_bc.evaluate(variation_str, working_point, f_sub[bc_idx], a_sub[bc_idx], p_sub[bc_idx])
            
            if np.any(l_idx):
                raw_sfs[l_idx] = corr_light.evaluate(variation_str, working_point, f_sub[l_idx], a_sub[l_idx], p_sub[l_idx])
                
            # Place the evaluated SFs back into the full array
            full_sfs[mask] = raw_sfs

        # ==========================================
        # --- PURE METHOD 1a MATH ---
        # ==========================================
        # Tagged jets get their true data SF
        weight_tagged = np.where(is_tagged, full_sfs, 1.0)
        
        # 1. Floating-Point Safety Net
        # We clip to 0.99999 strictly to prevent Numpy from throwing a Divide-by-Zero 
        # exception in statistically starved bins that hit exactly 1.0.
        safe_eff = np.clip(eff, 0.0, 0.99999)
        
        # 2. Pure Failing Weight Calculation (No Jet-Level Clamps)
        weight_untagged = np.where(~is_tagged, (1.0 - (full_sfs * safe_eff)) / (1.0 - safe_eff), 1.0)
        
        # 3. Combine passing and failing jets
        jet_weights = weight_tagged * weight_untagged
        
        # 4. Total Event Weight (No Event-Level Clamps)
        event_weights = ak.prod(ak.unflatten(jet_weights, njets), axis=1)
        
        return event_weights

    output = {}
    all_jets_mask = np.ones_like(flavour, dtype=bool)

    for variation in variations:
        if variation == "central":
            output[variation] = [_get_method1a_weight("central", all_jets_mask)]
        else:
            nominal = np.ones(ak.num(njets, axis=0)) 
            
            if variation in ["up", "down"]:
                output[variation] = [_get_method1a_weight(variation, all_jets_mask)]
                
            elif "cferr" in variation:
                c_mask = flavour == 4
                output[variation] = [
                    nominal,
                    _get_method1a_weight(f"up_{variation}", c_mask),
                    _get_method1a_weight(f"down_{variation}", c_mask),
                ]

            elif variation.startswith("JES") and "AK4" in variation:
                btag_jes_var = variation.replace("_AK4PFchs", "").replace("_AK4PFPuppi", "")
                if btag_jes_var.startswith("JES_Total") and btag_jes_var.endswith("Up"):
                    btag_jes_var = "up_jes"
                elif btag_jes_var.startswith("JES_Total") and btag_jes_var.endswith("Down"):
                    btag_jes_var = "down_jes"
                else:
                    if btag_jes_var.endswith("Up"):
                        btag_jes_var = f"up_jes{btag_jes_var[4:-2]}"
                    elif btag_jes_var.endswith("Down"):
                        btag_jes_var = f"down_jes{btag_jes_var[4:-4]}"

                # Evaluate JES variations only on non-c jets
                notc_mask = flavour != 4 
                output[variation] = [_get_method1a_weight(btag_jes_var, notc_mask)]

    return output
    
def apply_btag_sf(params, events, year, eff_lookup, working_point="M"):
    """
    Wrapper function to calculate b-tagging Method 1a scale factors.
    Handles data safety, flavor cleaning, and efficiency map evaluation.
    """
    jets = events.JetGood
    
    # --- 1. DATA SAFETY CHECK ---
    # If this is Data, return neutral event weights (1.0)
    if "hadronFlavour" not in jets.fields:
        ones = np.ones(len(events))
        return ones, ones, ones

    # --- 2. PREPARE KINEMATICS ---
    # Clean the flavor array so it strictly matches [0, 4, 5]
    raw_flavor = jets.hadronFlavour
    clean_flavour = ak.where((raw_flavor != 5) & (raw_flavor != 4), 0, raw_flavor)
    
    # Evaluate the lookup tool to get MC efficiencies
    mc_efficiencies = eff_lookup(
        jets.pt, 
        np.abs(jets.eta), 
        clean_flavour
    )
    
    # --- 3. APPLY SCALE FACTORS ---
    # Get the Working Point threshold from parameters
    wp_threshold = params.btagging.working_point[year]["btagging_WP"][working_point]
    
    # Call the core Method 1a math function
    btag_sf_dict = sf_btag_wp(
        params, 
        jets, 
        year, 
        njets=events.nJetGood, 
        mc_efficiencies=mc_efficiencies, 
        wp_threshold=wp_threshold,       
        variations=['central', 'up', 'down'], 
        working_point=working_point
    )
    
    return btag_sf_dict['central'][0], btag_sf_dict['up'][0], btag_sf_dict['down'][0]