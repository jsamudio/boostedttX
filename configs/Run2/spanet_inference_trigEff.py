import awkward as ak
import numpy as np
import vector
import correctionlib
import os

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
#from pocket_coffea.lib.weights_manager import WeightsManager
from object_cleaning_functions import soft_lep_sel, lep_sel, bjet_sel, qjet_sel, lep_softlep_combo, jet_selection_custom, lep_selTrig_Mu
from custom_cut_functions import sortbyscore
from analysis_years import xbb_wp_cuts
from cand_helper import zh_helper
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig, match_tt_products, match_gen_bb_to_zh
#import dnn_model
#from applyDNN import applyDNN
from weight_handler import calc_weight, add_weights_to_ttbb
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.scale_factors import sf_ele_reco, sf_ele_id, sf_btag, sf_pileup_reweight, sf_L1prefiring
from dask.distributed import get_worker
from custom_scale_factors import sf_mu, sf_partonshower_isr, sf_partonshower_fsr, sf_scaleweights, apply_jms_jmr, sf_ele_trig
import cachetools
from get_controlvars import get_controlvars

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

vector.register_awkward()

class ZHbbSpanetProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        self.jec_pt_variation = self.workflow_options["jec_pt_variation"]
        self.jer_variation = self.workflow_options["jer_variation"]
        self.pnet_mass_cset = correctionlib.CorrectionSet.from_file(self.workflow_options["pnet_mass_json"])
        # Check for the new spanet_models dictionary
        if not "spanet_models" in self.workflow_options:
            raise ValueError("Key `spanet_models` not found in workflow options. Please specify a dictionary of ONNX models.")


    def apply_object_preselection(self, variation):
        # 1. List all triggers that MIGHT be missing in early eras
        expected_triggers = [
            'Ele32_WPTight_Gsf', 'Photon200', 'Ele115_CaloIdVT_GsfTrkIdT', 
            'Ele50_CaloIdVT_GsfTrkIdT_PFJet165', 'IsoMu24', 'Mu50', 
            'OldMu100', 'TkMu100'
        ]
        
        # 2. Check the HLT record. If a trigger is missing, safely inject a False array
        for trig in expected_triggers:
            if trig not in self.events.HLT.fields:
                self.events = ak.with_field(
                    self.events, 
                    ak.zeros_like(self.events.run, dtype=bool), 
                    where=["HLT", trig]
                )
        #soft e, e, soft mu, mu, jet, fatjet
        
        # Map pocket_coffea's variation string to the correctionlib syst input
        syst_jms = "nominal"
        syst_jmr = "nominal"
        
        if variation == "JMSUp":
            syst_jms = "up"
        elif variation == "JMSDown":
            syst_jms = "down"
        elif variation == "JMRUp":
            syst_jmr = "up"
        elif variation == "JMRDown":
            syst_jmr = "down"
            
        if self._isMC:
            # Map PocketCoffea year strings to correctionlib era keys
            year_map_pnet = {
                "2016_PreVFP": "2016APV",
                "2016_PostVFP": "2016",
                "2017": "2017",
                "2018": "2018",
            }
            pnet_era = year_map_pnet.get(self._year, self._year)

            # Evaluate the SFs using the mapped era
            jms_sf = self.pnet_mass_cset["FatJet_PNet_JMS"].evaluate(pnet_era, syst_jms)
            jmr_sf = self.pnet_mass_cset["FatJet_PNet_JMR"].evaluate(pnet_era, syst_jmr)
            
            # The JMAR team utilized a 10% baseline for the MC resolution in their derivations
            sigma_mc = 0.10 
            
            # Apply corrections to the RAW FatJet collection before cuts
            corrected_mass = apply_jms_jmr(
                self.events, 
                self.events.FatJet, 
                jms_sf, 
                jmr_sf, 
                sigma_mc
            )
            
            self.events['FatJet'] = ak.with_field(
                self.events.FatJet, 
                corrected_mass, 
                "mass"
            )
            
            
        # Just needed for plotting
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )

        self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        #self.events['MuonGood'] = lep_selTrig_Mu(self.events, "Muon", self.params)
        self.events['SoftMuonGood'] = soft_lep_sel(self.events, "Muon", self.params)
        self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)
        self.events['SoftElectronGood'] = soft_lep_sel(self.events, "Electron", self.params)

        leptons = ak.with_name(
                ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis = 1),
                name='PtEtaPhiMCandidate')
        self.events['LeptonGood'] = leptons[ak.argsort(leptons.pt, ascending=False)]
        self.events['JetGood'], self.jetGoodMask = jet_selection_custom(self.events, "Jet", self.params, self._year, "LeptonGood", self.params.object_preselection.btag.tagger) #use PC default; should be deepFlav
        self.events['FatJetGood'], self.fatJetGoodMask = jet_selection_custom(self.events, "FatJet", self.params, self._year, "LeptonGood", "PNetMD") # they use some MSD cut which we don't want
        # Override the default 'mass' field so the vector library uses the regressed mass natively
        self.events['bJetGood'] = bjet_sel(self.events.JetGood, self.params, self._year)
        self.events['qJetGood'] = qjet_sel(self.events.JetGood, self.params, self._year)
        self.events['e_softe'] = lep_softlep_combo(self.events.ElectronGood, self.events.SoftElectronGood)
        self.events['mu_softmu'] = lep_softlep_combo(self.events.MuonGood, self.events.SoftMuonGood)

        # add columns for SPANet
        self.events['JetGood_pt'] = self.events.JetGood.pt
        self.events['JetGood_eta'] = self.events.JetGood.eta
        self.events['JetGood_phi'] = self.events.JetGood.phi
        #self.events['JetGood_btagDeepFlavB'] = self.events.JetGood.btagDeepFlavB
        #FIXME I think the btag LMH are working point boolean flags? Maybe omit for now
        self.events['LeptonGood_pt'] = self.events.LeptonGood.pt
        self.events['LeptonGood_eta'] = self.events.LeptonGood.eta
        self.events['LeptonGood_phi'] = self.events.LeptonGood.phi
        #self.events['MET_pt'] = self.events.MET.pt
        #self.events['MET_eta'] = self.events.LeptonGood.eta
        #self.events['MET_phi'] = self.events.MET.phi

    def do_parton_matching(self) -> ak.Array:
        # Selects quarks at LHE level
        isOutgoing = self.events.LHEPart.status == 1
        isParton = (abs(self.events.LHEPart.pdgId) < 6) | (self.events.LHEPart.pdgId == 21)
        quarks = self.events.LHEPart[isOutgoing & isParton]
        #print("quarks", ak.num(quarks, axis=-1))

        # Select b-quarks at Gen level, coming from Z/H->bb decay
        # for now seeing if we reuse the higgs collection for both
        if self._sample in ['ttHTobb', 'TTZToBB']:
            higgs = self.events.GenPart[
                ((self.events.GenPart.pdgId == 25) | (self.events.GenPart.pdgId == 23))
                & (self.events.GenPart.hasFlags(['fromHardProcess']))
            ]
            higgs = higgs[ak.num(higgs.childrenIdxG, axis=2) == 2]
            self.events["HiggsGen"] = higgs

            matched_higgs, matched_fatjets, deltaR_matchedAK8 = object_matching(
                higgs, self.events.FatJetSorted, dr_min=self.dr_min
            )
            #print("Matched fatjets:", matched_fatjets)

            if self._sample in ['ttHTobb']:
                higgs_partons = ak.with_field(
                    ak.flatten(higgs.children, axis=2), 25, "from_part"
                )
            else:
                higgs_partons = ak.with_field(
                    ak.flatten(higgs.children, axis=2), 23, "from_part"
                )
            # DO NOT sort b-quarks by pt
            # if not we are not able to match them with the provenance
            quarks = ak.with_name(
                ak.concatenate((quarks, higgs_partons), axis=1),
                name='PtEtaPhiMCandidate',
            )
        else:
            higgs = self.events.GenPart[
                (self.events.GenPart.hasFlags(['fromHardProcess']))
            ]
        #print("quarks", quarks.pdgId)
        # Get the interpretation
        if self._sample in ['ttHTobb', 'ttHTobb_ttToSemiLep']:
            prov = get_partons_provenance_ttHbb(
                ak.Array(quarks.pdgId, behavior={}), ak.ArrayBuilder()
            ).snapshot()
            self.events["HiggsParton"] = self.events.LHEPart[
                self.events.LHEPart.pdgId == 25
            ]
            higgs = higgs[higgs.pdgId == 25]
            ak8prov = 1 * ak.ones_like(higgs.pt)
        elif self._sample in ['TTZToBB']:
            prov = get_partons_provenance_ttHbb(
                ak.Array(quarks.pdgId, behavior={}), ak.ArrayBuilder()
            ).snapshot()
            self.events["HiggsParton"] = self.events.LHEPart[
                self.events.LHEPart.pdgId == 23
            ]
            higgs = higgs[higgs.pdgId == 23]
            ak8prov = 1 * ak.ones_like(higgs.pt)
        elif self._sample == "TTbb_SemiLeptonic":
            prov = get_partons_provenance_ttbb4F(
                ak.Array(quarks.pdgId, behavior={}), ak.ArrayBuilder()
            ).snapshot()
            ak8prov = -1 * ak.ones_like(higgs.pt)
        elif self._sample == "TTToSemiLeptonic":
            prov = get_partons_provenance_tt5F(
                ak.Array(quarks.pdgId, behavior={}), ak.ArrayBuilder()
            ).snapshot()
            ak8prov = -1 * ak.ones_like(higgs.pt)
        else:
            prov = -1 * ak.ones_like(quarks)
            ak8prov = -1 * ak.ones_like(higgs.pt)

      
        # Adding the provenance to the quark object
        quarks = ak.with_field(quarks, prov, "provenance")
        higgs = ak.with_field(higgs, ak8prov, "provenance")
        #print("Higgs:", higgs.provenance)
        #self.events["HiggsMatched"] = higgs

        # Calling our general object_matching function.
        # The output is an awkward array with the shape of the second argument and None where there is no matching.
        # So, calling like this, we will get out an array of matched_quarks with the dimension of the JetGood.
        matched_quarks, matched_jets, deltaR_matched = object_matching(
            quarks, self.events.JetGood, dr_min=self.dr_min
        )     
        matched_higgs, matched_fatjets, deltaR_matchedAK8 = object_matching(
            higgs, self.events.FatJetSorted, dr_min=0.6
        )
        self.events["FatJetMatched"] = ak.with_field(
            matched_fatjets, deltaR_matchedAK8, "dRMatchedJet"
        )
        #print("Higgs:", matched_higgs)
        self.events["HiggsMatched"] = matched_higgs
                              
        # Saving leptons and neutrino parton level
        self.events["LeptonParton"] = self.events.LHEPart[
            (self.events.LHEPart.status == 1)
            & (abs(self.events.LHEPart.pdgId) > 10)
            & (abs(self.events.LHEPart.pdgId) < 15)
        ]

        self.events["Parton"] = quarks
        self.events["PartonMatched"] = ak.with_field(
            matched_quarks, deltaR_matched, "dRMatchedJet"
        )
        self.events["JetGoodMatched"] = ak.with_field(
            matched_jets, deltaR_matched, "dRMatchedJet"
        )
        self.matched_partons_mask = ~ak.is_none(self.events.JetGoodMatched, axis=1)

    def count_partons(self):
        self.events["nParton"] = ak.num(self.events.Parton, axis=1)
        self.events["nPartonMatched"] = ak.count(
            self.events.PartonMatched.pt, axis=1
        )  # use count since we have None
        
    def get_ttbb_LHE_info(self):
        # Select the top quarks at LHE level
        isOutgoing = self.events.LHEPart.status == 1
        self.events['LHE_HT'] = self.events.LHE.HT
        # Events look like: [5, -5, 21, 5, -5, -3, 4, 15, -16] where the first two are the g->bb
        # Third particle is additional radiation
        # first top here is top, second is antitop
        # the first top gets the first string after the b quarks, so -s +c
        # thus it decayed hadronically
        # antitop then is leptonically decaying
        # so for this kind of event:
        # if len(pdgIds) == 10:
        #   bb = LHEPart[0] + LHEPart[1]
        #   t1 = LHEPart[3] + LHEPart[5] + LHEPart[6]
        #   t2 = LHEPart[4]+ LHEPart[7] + LHEPart[8]
        print("Outgoing: ", self.events.LHEPart[isOutgoing].pdgId)
        bb_ht = abs(self.events.LHEPart[isOutgoing].pt[:,0]) + abs(self.events.LHEPart[isOutgoing].pt[:,1])
        tt_ht = (abs(self.events.LHEPart[isOutgoing].pt[:,3]) + abs(self.events.LHEPart[isOutgoing].pt[:,4]) +
                 abs(self.events.LHEPart[isOutgoing].pt[:,5]) + abs(self.events.LHEPart[isOutgoing].pt[:,6]))# +
                 #abs(self.events.LHEPart[isOutgoing].pt[:,7]))# + abs(self.events.LHEPart[isOutgoing].pt[:,8]))
        glu_ht = abs(self.events.LHEPart[isOutgoing].pt[:,2])
        tt_withlep_ht = (abs(self.events.LHEPart[isOutgoing].pt[:,3]) + abs(self.events.LHEPart[isOutgoing].pt[:,4]) +
                 abs(self.events.LHEPart[isOutgoing].pt[:,5]) + abs(self.events.LHEPart[isOutgoing].pt[:,6]) +
                 abs(self.events.LHEPart[isOutgoing].pt[:,7]) + abs(self.events.LHEPart[isOutgoing].pt[:,8]))
        ttbb_full_ht = tt_withlep_ht + bb_ht
        ttbb_had_ht = tt_ht + bb_ht
        self.events["ttbb_full_ht"] = ttbb_full_ht
        self.events["ttbb_had_ht"] = ttbb_had_ht
        self.events["tt_ht"] = tt_ht
        self.events["bb_ht"] = bb_ht
        print("bbht", bb_ht)
        print("ttht", tt_ht)
        print("gluht", glu_ht)

    def define_common_variables_after_presel(self, variation):
        super().define_common_variables_before_presel(variation=variation)

        # Define labels for btagged jets at different working points
        for wp, val in self.params.btagging.working_point[self._year]["btagging_WP"].items():
            self.events["JetGood"] = ak.with_field(
                self.events.JetGood,
                ak.values_astype(self.events.JetGood.btagDeepFlavB > val, int),
                f"btag_{wp}"
            )
                # Compute the `is_electron` flag for LeptonGood
        self.events["LeptonGood"] = ak.with_field(
            self.events.LeptonGood,
            ak.values_astype(self.events.LeptonGood.pdgId == 11, bool),
            "is_electron"
        )
        xbbVsQCD = self.events.FatJetGood.particleNetMD_Xbb / (self.events.FatJetGood.particleNetMD_Xbb + self.events.FatJetGood.particleNetMD_QCD)
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, xbbVsQCD, 'xbbVsQCD')
        # --- NEW: Add boolean flags for the particleNetMD_XbbvsQCD Working Points ---
        for wp, cut in xbb_wp_cuts(self._year).items():   # analysis_years.XBB_WP_DICT
            self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, ak.values_astype(xbbVsQCD > cut, int), f'xbb_wp_{wp}')
        #self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, '2017')
        
    def calculate_xbb_sfs(self):
        import correctionlib
        import pickle
        import numpy as np
        
        # 1. Lazy-load evaluators and efficiency maps ONCE per worker
        if not hasattr(self, 'xbb_real_eval'):
            year_map = {"2016_PreVFP": "2016_PreVFP", "2016_PostVFP": "2016", "2017": "2017", "2018": "2018"}
            era = year_map.get(self._year, self._year)
            
            # Update these paths to your actual file locations
            real_sf_path = f"/cms/data/jsamudio/boosted/boostedttX/configs/Run2/SFBDT_Xbb_{era}.json"
            fake_sf_path = f"/cms/data/jsamudio/boosted/boostedttX/configs/Run2/Fake_SF_Xbb_{era}.json"
            eff_map_path = f"/cms/data/jsamudio/boosted/boostedttX/configs/Run2/efficiency_maps/efficiency_maps_{era}.pkl"
            
            self.xbb_real_eval = correctionlib.CorrectionSet.from_file(real_sf_path)[f"SFBDT_Xbb_SF_{era}"]
            self.xbb_fake_eval = correctionlib.CorrectionSet.from_file(fake_sf_path)[f"Fake_SF_Xbb_{era}"]
            with open(eff_map_path, 'rb') as f:
                self.eff_maps = pickle.load(f)

        # 2. Extract kinematics as flat NumPy arrays for fast processing
        fatjets = self.events.FatJetSorted[:, 0]
        pt = ak.to_numpy(fatjets.pt)
        msd = ak.to_numpy(fatjets.mass)
        score = ak.to_numpy(fatjets.xbbVsQCD)
        genmatch = ak.to_numpy(self.events.xbb_genmatch)
        n_events = len(pt)
 
        # 3. Discretize into mutually exclusive categories
        THRESH_HP = self.eff_maps['thresholds']['HP']
        THRESH_MP = self.eff_maps['thresholds']['MP']
        THRESH_LP = self.eff_maps['thresholds']['LP']
 
        is_hp = score >= THRESH_HP
        is_mp = (score >= THRESH_MP) & (score < THRESH_HP)
        is_lp = (score >= THRESH_LP) & (score < THRESH_MP)
        is_fail = score < THRESH_LP
 
        # Map pT and mSD to the efficiency dictionary indices
        pt_idx = np.clip(np.digitize(pt, self.eff_maps['pt_bins']) - 1, 0, len(self.eff_maps['pt_bins']) - 2)
        msd_idx = np.clip(np.digitize(msd, self.eff_maps['msd_bins']) - 1, 0, len(self.eff_maps['msd_bins']) - 2)
 
        # Process label PER EVENT, not per chunk -- match_gen_sig() can
        # assign different labels to different events within the same
        # sample (e.g. TTZToQQ splits into 'ttZ' vs a legacy 'old_ttZbb'
        # category for some events, presumably from an old overlap-removal
        # scheme against the dedicated TTZToBB sample). A single
        # representative label for the whole chunk would silently apply
        # the wrong per-process efficiency to whichever events don't
        # match it.
        #
        # Legacy label aliasing: 'old_ttZbb' events are physically
        # ordinary ttZ real-bb candidates (same tt-Z production, same
        # Z->bb decay) -- there's no physical reason their tagging
        # efficiency should differ from the rest of 'ttZ'. Whatever
        # overlap-removal/discard step already exists downstream still
        # removes these events from the final accounting regardless of
        # what weight they carry here, since the weight and the discard
        # are independent per-event operations. Extend this dict if
        # other legacy labels turn up.
        PROCESS_ALIASES = {
            'old_ttZbb': 'ttZ',
        }
        process_arr = ak.to_numpy(self.events['process']).copy()
        for old_label, new_label in PROCESS_ALIASES.items():
            process_arr[process_arr == old_label] = new_label
 
        per_process_effmaps = self.eff_maps['real']['per_process']
        known_keys = set(per_process_effmaps.keys())
 
        # Safety net for any label that's genuinely unrecognized even
        # after aliasing -- still defaults to an inert SF=1.0 with a
        # visible warning, rather than crashing the whole file.
        unknown_mask = ~np.isin(process_arr, list(known_keys))
        if np.any(unknown_mask):
            unknown_labels = np.unique(process_arr[unknown_mask])
            print(f"\u26A0\uFE0F  {int(unknown_mask.sum())} event(s) have process label(s) "
                  f"{list(unknown_labels)} not found in eff_maps['real']['per_process'] "
                  f"(known: {sorted(known_keys)}) even after PROCESS_ALIASES. Defaulting "
                  f"xbb_sf=1.0 for these events -- investigate whether this is a new "
                  f"legacy label that needs an alias, or a genuine labeling bug.")
 
        eff_hp = np.ones(n_events, dtype=np.float64)
        eff_mp = np.ones(n_events, dtype=np.float64)
        eff_lp = np.ones(n_events, dtype=np.float64)
        for proc in known_keys:
            proc_mask = (process_arr == proc)
            if not np.any(proc_mask):
                continue
            eff_hp[proc_mask] = per_process_effmaps[proc]["HP"]['eff'][pt_idx[proc_mask], msd_idx[proc_mask]]
            eff_mp[proc_mask] = per_process_effmaps[proc]["MP"]['eff'][pt_idx[proc_mask], msd_idx[proc_mask]]
            eff_lp[proc_mask] = per_process_effmaps[proc]["LP"]['eff'][pt_idx[proc_mask], msd_idx[proc_mask]]
 
        sf_results = {}
 
        # 4. Loop over systematic variations
        for syst in ["nominal", "up", "down"]:
            # --- REAL CANDIDATES (Method 1a Multi-WP) ---
            sf_real_hp = np.array([self.xbb_real_eval.evaluate("HP", syst, p) for p in pt])
            sf_real_mp = np.array([self.xbb_real_eval.evaluate("MP", syst, p) for p in pt])
            sf_real_lp = np.array([self.xbb_real_eval.evaluate("LP", syst, p) for p in pt])
 
            w_real_hp = sf_real_hp
 
            num_mp = (sf_real_mp * eff_mp) - (sf_real_hp * eff_hp)
            den_mp = eff_mp - eff_hp
            w_real_mp = np.where(num_mp < 0, 1.0, num_mp / np.where(den_mp > 0, den_mp, 1.0))
 
            num_lp = (sf_real_lp * eff_lp) - (sf_real_mp * eff_mp)
            den_lp = eff_lp - eff_mp
            w_real_lp = np.where(num_lp < 0, 1.0, num_lp / np.where(den_lp > 0, den_lp, 1.0))
 
            real_weights = np.ones_like(pt)
            real_weights[is_hp] = w_real_hp[is_hp]
            real_weights[is_mp] = w_real_mp[is_mp]
            real_weights[is_lp] = w_real_lp[is_lp]
 
            # --- FAKE CANDIDATES (Data-Driven) ---
            fake_weights = np.ones_like(pt)
            if np.any(is_hp):   fake_weights[is_hp]   = np.array([self.xbb_fake_eval.evaluate("HP", syst, p, m) for p, m in zip(pt[is_hp], msd[is_hp])])
            if np.any(is_mp):   fake_weights[is_mp]   = np.array([self.xbb_fake_eval.evaluate("MP", syst, p, m) for p, m in zip(pt[is_mp], msd[is_mp])])
            if np.any(is_lp):   fake_weights[is_lp]   = np.array([self.xbb_fake_eval.evaluate("LP", syst, p, m) for p, m in zip(pt[is_lp], msd[is_lp])])
            if np.any(is_fail): fake_weights[is_fail] = np.array([self.xbb_fake_eval.evaluate("Fail", syst, p, m) for p, m in zip(pt[is_fail], msd[is_fail])])
 
            # --- SPLICE TOGETHER ---
            final_sf = np.ones_like(pt)
            final_sf[genmatch] = real_weights[genmatch]
            final_sf[~genmatch] = fake_weights[~genmatch]
 
            # Hard override: events with an unrecognized process label
            # (e.g. 'old_ttZbb') always get an inert 1.0, regardless of
            # whatever real_weights/fake_weights computed above using
            # the default eff=1.0 fallback -- that fallback keeps eff_mp
            # den_mp etc. from producing garbage, but doesn't by itself
            # guarantee a clean no-op weight through the formula.
            final_sf[unknown_mask] = 1.0
 
            sf_results[syst] = ak.from_numpy(final_sf)
 
        return sf_results["nominal"], sf_results["up"], sf_results["down"]



    def process_extra_after_presel(self, variation):
        #FIXME
        
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")
        vec4 = vector.zip(
            {
                "pt": self.events.FatJetSorted[:,0].pt,
                "eta": self.events.FatJetSorted[:,0].eta,
                "phi": self.events.FatJetSorted[:,0].phi,
                "mass": self.events.FatJetSorted[:,0].mass,
            }
        )
        self.events['ZH_rap'] = vec4.rapidity
        print(self.events.ZH_rap)
        
        # --- MC ONLY SECTION ---
        if 'DATA' not in self._sample:
            self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, self._year)
            print("ele_reco", self.events.ele_reco_sf)
            self.events['ele_id_sf'], self.events['ele_id_sfup'], self.events['ele_id_sfdown'] = sf_ele_id(self.params, self.events, self._year)
            print("ele_id", self.events.ele_id_sf)
            print("ele eta", self.events.ElectronGood.eta)
            self.events['ele_trig_sf'], self.events['ele_trig_sfup'], self.events['ele_trig_sfdown'] = sf_ele_trig(self.events, f"/cms/data/jsamudio/boosted/boostedttX/configs/Run2/ele_trig_sf_{self._year}.json.gz")
            print("ele_trig", self.events.ele_trig_sf)
            self.events['mu_id_sf'], self.events['mu_id_sfup'], self.events['mu_id_sfdown'], self.events['mu_id_sfstat'] = sf_mu(self.params, self.events, self._year, 'id')
            print("mu_id", self.events.mu_id_sf)
            self.events['mu_iso_sf'], self.events['mu_iso_sfup'], self.events['mu_iso_sfdown'], self.events['mu_iso_sfstat'] = sf_mu(self.params, self.events, self._year, 'iso')
            print("mu_iso", self.events.mu_iso_sf)
            self.events['mu_reco_sf'], self.events['mu_reco_sfup'], self.events['mu_reco_sfdown'], self.events['mu_reco_sfstat'] = sf_mu(self.params, self.events, self._year, 'reco')
            print("mu_reco", self.events.mu_reco_sf)
            self.events['mu_trig_sf'], self.events['mu_trig_sfup'], self.events['mu_trig_sfdown'], self.events['mu_trig_sfstat'] = sf_mu(self.params, self.events, self._year, 'trigger')
            print("mu_reco", self.events.mu_trig_sf)
            
            btag_variations = self.params.systematic_variations.weight_variations.sf_btag[self._year]
            btag_sf = sf_btag(self.params, self.events.JetGood, self._year, njets=self.events.nJetGood, variations=['central']+btag_variations)
            
            # btagging
            self.events['btag_sf'] = btag_sf['central'][0]
            self.events['btag_sfhf'] = btag_sf['hf'][0]
            self.events['btag_sfhf_up'] = btag_sf['hf'][1]
            self.events['btag_sfhf_down'] = btag_sf['hf'][2]
            self.events['btag_sflf'] = btag_sf['lf'][0]
            self.events['btag_sflf_up'] = btag_sf['lf'][1]
            self.events['btag_sflf_down'] = btag_sf['lf'][2]
            self.events['btag_sfhfstats1'] = btag_sf['hfstats1'][0]
            self.events['btag_sfhfstats1_up'] = btag_sf['hfstats1'][1]
            self.events['btag_sfhfstats1_down'] = btag_sf['hfstats1'][2]
            self.events['btag_sfhfstats2'] = btag_sf['hfstats2'][0]
            self.events['btag_sfhfstats2_up'] = btag_sf['hfstats2'][1]
            self.events['btag_sfhfstats2_down'] = btag_sf['hfstats2'][2]
            self.events['btag_sflfstats1'] = btag_sf['lfstats1'][0]
            self.events['btag_sflfstats1_up'] = btag_sf['lfstats1'][1]
            self.events['btag_sflfstats1_down'] = btag_sf['lfstats1'][2]
            self.events['btag_sflfstats2'] = btag_sf['lfstats2'][0]
            self.events['btag_sflfstats2_up'] = btag_sf['lfstats2'][1]
            self.events['btag_sflfstats2_down'] = btag_sf['lfstats2'][2]
            self.events['btag_sfcferr1'] = btag_sf['cferr1'][0]
            self.events['btag_sfcferr1_up'] = btag_sf['cferr1'][1]
            self.events['btag_sfcferr1_down'] = btag_sf['cferr1'][2]
            self.events['btag_sfcferr2'] = btag_sf['cferr2'][0]
            self.events['btag_sfcferr2_up'] = btag_sf['cferr2'][1]
            self.events['btag_sfcferr2_down'] = btag_sf['cferr2'][2]
            
            print("BTAG SF: ", len(btag_sf['central'][0]))
            print("len of events", len(self.events.ele_id_sf))
            print("number of jets", self.events.nJetGood)
            
            self.events['puWeight'], self.events['puWeight_up'], self.events['puWeight_down'] = sf_pileup_reweight(self.params, self.events, self._year)
            self.events['isr'], self.events['isr_up'], self.events['isr_down'] = sf_partonshower_isr(self.events)
            self.events['fsr'], self.events['fsr_up'], self.events['fsr_down'] = sf_partonshower_fsr(self.events)
            self.events['mu_r_up'], self.events['mu_r_down'], self.events['mu_f_up'], self.events['mu_f_down'], self.events['mu_rf_up'], self.events['mu_rf_down'], = sf_scaleweights(self.events)
            self.events['l1prefire_sf'], self.events['l1prefire_sfup'], self.events['l1prefire_sfdown'] = sf_L1prefiring(self.events)
            
            print("fsr_up", self.events.fsr_up)
            print("isr_up", self.events.isr_up)
            
            #match_gen_lep(self.events)

            match_gen_bb_to_zh(self.events)
            
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        elif (("TTbb" in self._sample) | ("TTTo" in self._sample)):
            match_gen_tt(self.events, self._sample)
        elif "JetsToL" in self._sample:
            self.events['process'] = 'VJets'
            self.events['topptWeight'] = 1
            self.events['topptWeight_Up'] = 1
            self.events['topptWeight_Down'] = 1
        else:
            self.events['process'] = self._sample
            self.events['topptWeight'] = 1
            self.events['topptWeight_Up'] = 1
            self.events['topptWeight_Down'] = 1
            
        if 'TTbb' in self._sample:
            add_weights_to_ttbb(self.events, self._sample)



        # --- APPLIES TO BOTH MC AND DATA ---
        
        #zh_helper(self.events)
        
        calc_weight(self.events, self.output, self._dataset, self.params, self._year)
        
        get_controlvars(self.events)

        #print(self.events.HLT.fields)

        #self.onnx_inference()

        # ==========================================================
        # XBB TAGGER SCALE FACTORS
        # ==========================================================
        # if ('DATA' not in self._sample):
        #     xbb_nom, xbb_up, xbb_dn = self.calculate_xbb_sfs()
        #     self.events['xbb_sf'] = xbb_nom
        #     self.events['xbb_sfup'] = xbb_up
        #     self.events['xbb_sfdown'] = xbb_dn
        # # ==========================================================

        # print("Double-B SF:", self.events.xbb_sf)

    def count_objects(self, variation):
        self.events['nMuonGood'] = ak.num(self.events.MuonGood)
        self.events['nElectronGood'] = ak.num(self.events.ElectronGood)
        self.events['nSoftMuonGood'] = ak.num(self.events.SoftMuonGood)
        self.events['nSoftElectronGood'] = ak.num(self.events.SoftElectronGood)
        self.events['nJetGood'] = ak.num(self.events.JetGood)
        self.events['nFatJetGood'] = ak.num(self.events.FatJetGood)
        self.events['nFatJetGood2'] = ak.num(self.events.FatJetGood)
        self.events['nbJetGood'] = ak.num(self.events.bJetGood)
        self.events['nLeptonGood'] = (self.events['nMuonGood'] + self.events['nElectronGood'])

    def onnx_inference(self) -> ak.Array:
        import onnxruntime as ort
        
        try:
            worker = get_worker()
        except ValueError:
            worker = None

        # ==========================================================
        # 1. PREPARE TWO JET COLLECTIONS: FULL & CLEANED
        # ==========================================================
        # Full collection
        jets_full = self.events.JetGood
        
        # Cleaned collection (Overlap Removal)
        leading_fatjet = ak.firsts(self.events.FatJetSorted)
        dr_jets = self.events.JetGood.delta_r(leading_fatjet)
        overlap_mask = (dr_jets > 0.8) | ak.is_none(dr_jets)
        jets_cleaned = self.events.JetGood[overlap_mask]

        # Helper function to pad and format jet variables
        def get_jet_features(jet_collection):
            pad_dict = {"btagB":0., "pt":0., "phi":0., "eta":0.}
            jets_padded = ak.zip(
                {key : ak.fill_none(ak.pad_none(jet_collection[key], 16, clip=True), value) for key, value in pad_dict.items()}
            )
            jet_vars = {
                "pt": np.log(1 + ak.to_numpy(jets_padded.pt)),
                "eta": ak.to_numpy(jets_padded.eta),
                "sin_phi": np.sin(ak.to_numpy(jets_padded.phi)),
                "cos_phi": np.cos(ak.to_numpy(jets_padded.phi)),
                "btagB": ak.to_numpy(jets_padded.btagB)
            }
            mask = ~ak.to_numpy(jets_padded.pt == 0)
            return jet_vars, mask

        # Extract numpy arrays for both sets
        jet_vars_full, mask_full = get_jet_features(jets_full)
        jet_vars_cleaned, mask_cleaned = get_jet_features(jets_cleaned)

        # ==========================================================
        # 2. PREPARE STATIC FEATURES (MET, Lepton, FatJet)
        # ==========================================================
        lep_vars = {
            "pt": np.log(1 + ak.to_numpy(self.events.LeptonGood[:,0].pt)),
            "eta": ak.to_numpy(self.events.LeptonGood[:,0].eta),
            "sin_phi": np.sin(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
            "cos_phi": np.cos(ak.to_numpy(self.events.LeptonGood[:,0].phi))
        }
        
        met_vars = {
            "pt": np.log(1+ ak.to_numpy(self.events.MET.pt)),
            "sin_phi": np.sin(ak.to_numpy(self.events.MET.phi)),
            "cos_phi": np.cos(ak.to_numpy(self.events.MET.phi))
        }

        fatjet_vars = {
            "pt": np.log(1 + ak.to_numpy(self.events.FatJetSorted[:,0].pt)),
            "eta": ak.to_numpy(self.events.FatJetSorted[:,0].eta),
            "sin_phi": np.sin(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
            "cos_phi": np.cos(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
            "bbtag": ak.to_numpy(self.events.FatJetSorted[:,0].xbbVsQCD),
            "xbb_wp_L": ak.to_numpy(self.events.FatJetSorted[:,0].xbb_wp_L),
            "xbb_wp_M": ak.to_numpy(self.events.FatJetSorted[:,0].xbb_wp_M),
            "xbb_wp_T": ak.to_numpy(self.events.FatJetSorted[:,0].xbb_wp_T)
        }

        mask_global = np.ones(shape=[met_vars["pt"].shape[0], 1]) == 1

        # ==========================================================
        # 3. LOOP OVER DEFINED MODELS
        # ==========================================================
        for model_name, model_info in self.workflow_options["spanet_models"].items():
            model_file = model_info["path"]
            
            # --- Session Caching on Dask Workers ---
            if worker is not None:
                if getattr(worker, 'spanet_sessions', None) is None:
                    worker.spanet_sessions = {}
                if model_name not in worker.spanet_sessions:
                    sess_options = ort.SessionOptions()
                    sess_options.intra_op_num_threads = 1
                    sess_options.inter_op_num_threads = 1
                    sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
                    worker.spanet_sessions[model_name] = ort.InferenceSession(model_file, sess_options=sess_options, providers=['CPUExecutionProvider'])
                model_session = worker.spanet_sessions[model_name]
            else:
                sess_options = ort.SessionOptions()
                model_session = ort.InferenceSession(model_file, sess_options=sess_options, providers=['CPUExecutionProvider'])
            
            # --- Dynamically Select Jet Array Based on Model Name ---
            if "NoOverlap" in model_name:
                current_jet_vars = jet_vars_cleaned
                current_jet_mask = mask_cleaned
            else:
                current_jet_vars = jet_vars_full
                current_jet_mask = mask_full
            
            # --- Stack Requested Features ---
            req_jet = model_info.get("jet_features", ["pt", "eta", "sin_phi", "cos_phi", "btagB"])
            req_fatjet = model_info.get("fatjet_features", ["pt", "eta", "sin_phi", "cos_phi", "bbtag"])
            
            jet_data = np.transpose(np.stack([current_jet_vars[f] for f in req_jet]), axes=[1,2,0]).astype(np.float32)
            fatjet_data = np.stack([fatjet_vars[f] for f in req_fatjet], axis=1)[:,None,:].astype(np.float32)
            
            met_data = np.stack([met_vars["pt"], met_vars["sin_phi"], met_vars["cos_phi"]], axis=1)[:,None,:].astype(np.float32)
            lep_data = np.stack([lep_vars["pt"], lep_vars["eta"], lep_vars["sin_phi"], lep_vars["cos_phi"]], axis=1)[:,None,:].astype(np.float32)

            # Define outputs
            if "Zonly" in model_file:
                output_names = ["EVENT/ttzbb"]
            else:
                output_names = ["EVENT/signal", "encoder_event_vector"]
                
            # Run inference
            outputs = model_session.run(
                input_feed={
                    "Jet_data": jet_data,
                    "Jet_mask": current_jet_mask,
                    "Met_data": met_data,
                    "Met_mask": mask_global,
                    "Lepton_data": lep_data,
                    "Lepton_mask": mask_global,
                    "FatJet_data": fatjet_data,
                    "FatJet_mask": mask_global
                },
                output_names=output_names
            )
            
            outputs_zipped = dict(zip(output_names, outputs))
            keys = list(outputs_zipped.keys())
            
            # 4. SAVE OUTPUTS WITH MODEL IDENTIFIER
            if "Zonly" in model_file:
                self.events[f"spanet_outputZ_{model_name}"] = ak.zip(
                    {key.split("/")[-1]: ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()}
                )
            else:
                self.events[f"spanet_output_{model_name}"] = ak.zip(
                    {keys[0].split("/")[-1]: ak.from_numpy(outputs_zipped[keys[0]][:,1]) / (ak.from_numpy(outputs_zipped[keys[0]][:,0]) + ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2]) + ak.from_numpy(outputs_zipped[keys[0]][:,3]))}
                )
                self.events[f"sig_score_{model_name}"] = ak.from_numpy(outputs_zipped[keys[0]][:,1])
                self.events[f"ttbb_score_{model_name}"] = ak.from_numpy(outputs_zipped[keys[0]][:,2])
                self.events[f"ttlf_score_{model_name}"] = ak.from_numpy(outputs_zipped[keys[0]][:,3])
                
                for i in range(32):
                    self.events[f"encoderFeature{i+1}_{model_name}"] = ak.from_numpy(outputs_zipped[keys[1]][:,i])

    # def onnx_inference(self, model_file) -> ak.Array:

    #     try:
    #         worker = get_worker()
    #     except ValueError:
    #         worker = None
            
    #     if worker is None:
    #         import onnxruntime as ort
    #         sess_options = ort.SessionOptions()
    #         sess_options.intra_op_num_threads = 1
    #         sess_options.inter_op_num_threads = 1
    #         sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
    #         model_session = ort.InferenceSession(
    #             model_file,
    #             sess_options = sess_options,
    #             providers=['CPUExecutionProvider']
    #         )
    #     else:
    #         model_session = worker.data['model_session']

    #     #print(model_session)

    #     btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]
    #     pad_dict = {"btagB":0., "pt":0., "phi":0., "eta":0.}
    #     jets_padded = ak.zip(
    #         {key : ak.fill_none(ak.pad_none(self.events.JetGood[key], 16, clip=True), value) for key, value in pad_dict.items()}
    #     )

    #     data = np.transpose(
    #         np.stack([
    #             np.log(1 + ak.to_numpy(jets_padded.pt)),
    #             ak.to_numpy(jets_padded.eta),
    #             np.sin(ak.to_numpy(jets_padded.phi)),
    #             np.cos(ak.to_numpy(jets_padded.phi)),
    #             ak.to_numpy(jets_padded.btagB),
    #         ]),
    #         axes=[1,2,0]).astype(np.float32)

    #     mask = ~ak.to_numpy(jets_padded.pt == 0)

    #     met_data = np.stack([np.log(1+ ak.to_numpy(self.events.MET.pt)), # pT
    #                          np.sin(ak.to_numpy(self.events.MET.phi)), #sin(phi)
    #                          np.cos(ak.to_numpy(self.events.MET.phi)) #cos(phi)
    #                          ], axis=1)[:,None,:].astype(np.float32)

    #     lep_data = np.stack([np.log(1 + ak.to_numpy(self.events.LeptonGood[:,0].pt)),
    #                          ak.to_numpy(self.events.LeptonGood[:,0].eta),
    #                          np.sin(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
    #                          np.cos(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
    #                          ], axis=1)[:,None,:].astype(np.float32)

    #     fatjet_data = np.stack([np.log(1 + ak.to_numpy(self.events.FatJetSorted[:,0].pt)),
    #                          ak.to_numpy(self.events.FatJetSorted[:,0].eta),
    #                          np.sin(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
    #                          np.cos(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
    #                          ak.to_numpy(self.events.FatJetSorted[:,0].btagBB),
    #                          ], axis=1)[:,None,:].astype(np.float32)
    #     event_data = ak.to_numpy(self.events.n_b_inZH[:,None,None]).astype(np.float32)

    #     mask_global = np.ones(shape=[met_data.shape[0], 1]) == 1
    #     if "Zonly" in model_file:
    #         output_names = ["EVENT/ttzbb"]
    #     elif "1" in model_file:
    #         output_names = ["EVENT/tthbb", "EVENT/ttbb", "EVENT/ttcc", "EVENT/ttlf"]
    #     #else:
    #     #output_names = ["EVENT/ttzbb", "EVENT/tthbb", "EVENT/ttbb", "EVENT/ttcc", "EVENT/ttlf"]
    #     #output_names = ["EVENT/ttzbb", "EVENT/tthbb", "EVENT/ttbb", "EVENT/ttlf"]
    #     output_names = ["EVENT/signal", "encoder_event_vector"]
    #     outputs = model_session.run(input_feed={
    #         "Jet_data": data,
    #         "Jet_mask": mask,
    #         "Met_data": met_data,
    #         "Met_mask": mask_global,
    #         "Lepton_data": lep_data,
    #         "Lepton_mask": mask_global,
    #         "FatJet_data": fatjet_data,
    #         "FatJet_mask": mask_global},
    #         #"Event_data": event_data,
    #         #"Event_mask": mask_global},
    #     output_names=output_names
    #     )
    #     #print("-" * 60)
    #     #print("OUTPUTS")
    #     #print(outputs)
    #     #print("-" * 60)
    #     #print((ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])))
       

    #     outputs_zipped = dict(zip(output_names, outputs))
    #     print(outputs_zipped.keys())
    #     keys = list(outputs_zipped.keys())
    #     #print((ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) for key, value in outputs_zipped.items())
    #     #print("0", (ak.from_numpy(value[:,0]) for key, value in outputs_zipped.items()))
    #     #print("1", (ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()))
    #     #print("2", (ak.from_numpy(value[:,2]) for key, value in outputs_zipped.items()))
    #     #print("3", (ak.from_numpy(value[:,3]) for key, value in outputs_zipped.items()))
    #     #print("4", (ak.from_numpy(value[:,4]) for key, value in outputs_zipped.items()))
    #     print(len(outputs_zipped[keys[1]][0,:]))
    #     if "Zonly" in model_file:
    #         self.events["spanet_outputZ"] = ak.zip(
    #             {
    #                 key.split("/")[-1]: ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()
    #             }
    #         )
    #     #else:
    #     self.events["spanet_output"] = ak.zip(
    #         {
    #             #key.split("/")[-1]:  (ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) for key, value in outputs_zipped.items()
    #             #keys[0].split("/")[-1]:  (ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2])) / (ak.from_numpy(outputs_zipped[keys[0]][:,0]) + ak.from_numpy(outputs_zipped[keys[0]][:,3]) + ak.from_numpy(outputs_zipped[keys[0]][:,4]) + ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2]))
    #             keys[0].split("/")[-1]:  ak.from_numpy(outputs_zipped[keys[0]][:,1]) / (ak.from_numpy(outputs_zipped[keys[0]][:,0]) + ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2]) + ak.from_numpy(outputs_zipped[keys[0]][:,3]))
    #             #"out_1": ak.from_numpy(value[:,1] for key, value in outputs_zipped.items()),
    #             #"out_2": ak.from_numpy(value[:,2] for key, value in outputs_zipped.items()),
    #         }
    #     )
    #     self.events["sig_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,1])
    #     self.events["ttbb_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,2])
    #     self.events["ttlf_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,3])
    #     self.events["encoderFeature1"] = ak.from_numpy(outputs_zipped[keys[1]][:,0])
    #     self.events["encoderFeature2"] = ak.from_numpy(outputs_zipped[keys[1]][:,1])
    #     self.events["encoderFeature3"] = ak.from_numpy(outputs_zipped[keys[1]][:,2])
    #     self.events["encoderFeature4"] = ak.from_numpy(outputs_zipped[keys[1]][:,3])
    #     self.events["encoderFeature5"] = ak.from_numpy(outputs_zipped[keys[1]][:,4])
    #     self.events["encoderFeature6"] = ak.from_numpy(outputs_zipped[keys[1]][:,5])
    #     self.events["encoderFeature7"] = ak.from_numpy(outputs_zipped[keys[1]][:,6])
    #     self.events["encoderFeature8"] = ak.from_numpy(outputs_zipped[keys[1]][:,7])
    #     self.events["encoderFeature9"] = ak.from_numpy(outputs_zipped[keys[1]][:,8])
    #     self.events["encoderFeature10"] = ak.from_numpy(outputs_zipped[keys[1]][:,9])
    #     self.events["encoderFeature11"] = ak.from_numpy(outputs_zipped[keys[1]][:,10])
    #     self.events["encoderFeature12"] = ak.from_numpy(outputs_zipped[keys[1]][:,11])
    #     self.events["encoderFeature13"] = ak.from_numpy(outputs_zipped[keys[1]][:,12])
    #     self.events["encoderFeature14"] = ak.from_numpy(outputs_zipped[keys[1]][:,13])
    #     self.events["encoderFeature15"] = ak.from_numpy(outputs_zipped[keys[1]][:,14])
    #     self.events["encoderFeature16"] = ak.from_numpy(outputs_zipped[keys[1]][:,15])
    #     self.events["encoderFeature17"] = ak.from_numpy(outputs_zipped[keys[1]][:,16])
    #     self.events["encoderFeature18"] = ak.from_numpy(outputs_zipped[keys[1]][:,17])
    #     self.events["encoderFeature19"] = ak.from_numpy(outputs_zipped[keys[1]][:,18])
    #     self.events["encoderFeature20"] = ak.from_numpy(outputs_zipped[keys[1]][:,19])
    #     self.events["encoderFeature21"] = ak.from_numpy(outputs_zipped[keys[1]][:,20])
    #     self.events["encoderFeature22"] = ak.from_numpy(outputs_zipped[keys[1]][:,21])
    #     self.events["encoderFeature23"] = ak.from_numpy(outputs_zipped[keys[1]][:,22])
    #     self.events["encoderFeature24"] = ak.from_numpy(outputs_zipped[keys[1]][:,23])
    #     self.events["encoderFeature25"] = ak.from_numpy(outputs_zipped[keys[1]][:,24])
    #     self.events["encoderFeature26"] = ak.from_numpy(outputs_zipped[keys[1]][:,25])
    #     self.events["encoderFeature27"] = ak.from_numpy(outputs_zipped[keys[1]][:,26])
    #     self.events["encoderFeature28"] = ak.from_numpy(outputs_zipped[keys[1]][:,27])
    #     self.events["encoderFeature29"] = ak.from_numpy(outputs_zipped[keys[1]][:,28])
    #     self.events["encoderFeature30"] = ak.from_numpy(outputs_zipped[keys[1]][:,29])
    #     self.events["encoderFeature31"] = ak.from_numpy(outputs_zipped[keys[1]][:,30])
    #     self.events["encoderFeature32"] = ak.from_numpy(outputs_zipped[keys[1]][:,31])
    #     print(self.events.encoderFeature20)
    #     #print(self.events.spanet_outputH)