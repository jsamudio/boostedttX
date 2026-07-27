import awkward as ak
import numpy as np
import vector

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
#from pocket_coffea.lib.weights_manager import WeightsManager
from pocket_coffea.lib.objects import (
    #jet_correction,
    lepton_selection,
    jet_selection,
    btagging,
    get_dilepton,
)
from object_cleaning_functions import soft_lep_sel, lep_sel, fatjet_sel, bjet_sel, qjet_sel, lep_softlep_combo, jet_sel, fatjet_sel2, jet_selection_custom
from custom_cut_functions import sortbyscore
from cand_helper import zh_helper
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig, match_tt_products
#import dnn_model
#from applyDNN import applyDNN
from weight_handler import calc_weight
from coffea.analysis_tools import PackedSelection
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.jets import compute_jetId, jet_selection
from pocket_coffea.lib.scale_factors import *
from dask.distributed import get_worker
from custom_scale_factors import sf_scaleweights, sf_btag_wp, sf_ele_trig, apply_btag_sf
from get_controlvars import get_controlvars
import cachetools

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu', 'ttHSMEFT', 'TTLL', 'TTNuNu']

vector.register_awkward()

class ZHbbSpanetProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        self.jec_pt_variation = self.workflow_options["jec_pt_variation"]
        self.jer_variation = self.workflow_options["jer_variation"]
        if not "spanet_model" in self.workflow_options:
            raise ValueError("Key `spanet_model` not found in workflow options. Please specify the path to the ONNX model.")
        elif not self.workflow_options["spanet_model"].endswith(".onnx"):
            raise ValueError("Key `spanet_model` should be the path of an ONNX model.")
        import pickle
        eff_map_path = "/cms/data/jsamudio/boosted/boostedttX/configs/SemiLeptonic_2024/BtagEff/semileptonic_btag_eff_2024.pkl"
        with open(eff_map_path, "rb") as f_eff:
            self.btag_eff_maps = pickle.load(f_eff)

    def skim_events(self):
        self._skim_masks = PackedSelection()
        mask_flags = np.ones(self.nEvents_initial, dtype=bool)

        flags = self.params.event_flags[self._year]
        if not self._isMC:
            flags += self.params.event_flags_data[self._year]
        for flag in flags:
            mask_flags &= getattr(self.events.Flag, flag).to_numpy()
        self._skim_masks.add("event_flags", mask_flags)

        for skim_func in self._skim:
            # Apply the skim function and add it to the mask
            mask = skim_func.get_mask(
                self.events,
                processor_params=self.params,
                year=self._year,
                sample=self._sample,
                isMC=self._isMC,
            )
            self._skim_masks.add(skim_func.id, mask)
        # Finally we skim the events and count them
        self.events = self.events[self._skim_masks.all(*self._skim_masks.names)]
        self.nEvents_after_skim = self.nevents
        self.output['cutflow']['skim'][self._dataset] = self.nEvents_after_skim
        self.has_events = self.nEvents_after_skim > 0

    def apply_object_preselection(self, variation):
        #soft e, e, soft mu, mu, jet, fatjet

        # Just needed for plotting
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(
            self.events.Electron, electron_etaSC, "etaSC"
        )

        self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        self.events['SoftMuonGood'] = soft_lep_sel(self.events, "Muon", self.params)
        self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)
        self.events['SoftElectronGood'] = soft_lep_sel(self.events, "Electron", self.params)

        leptons = ak.with_name(
                ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis = 1),
                name='PtEtaPhiMCandidate')
        self.events['LeptonGood'] = leptons[ak.argsort(leptons.pt, ascending=False)]

        self.events['JetGood'], self.jetGoodMask = jet_selection_custom(self.events, "Jet", self.params, self._year, "LeptonGood", self.params.object_preselection.btag.tagger) #use PC default
        self.events['FatJetGood'], self.fatJetGoodMask = jet_selection_custom(self.events, "FatJet", self.params, self._year, "LeptonGood", "GloParT") # they use some MSD cut which we don't want
        
        #print(self.events.JetGood.pt)
        #self.events.JetGood["hadronFlavour"] = ak.values_astype(self.events.JetGood["hadronFlavour"], "int32")
        
        #elf.events['FatJetGood'] = fatjet_sel(self.events, self.params, "LeptonGood")
        #self.events['FatJetGood2'] = fatjet_sel2(self.events, self.params, "LeptonGood")
        self.events['bJetGood'] = bjet_sel(self.events.JetGood, self.params)
        self.events['qJetGood'] = qjet_sel(self.events.JetGood, self.params)
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
        #self.events["LeptonGood"] = ak.with_field(
        #    self.events.LeptonGood,
        #    ak.values_astype(self.events.LeptonGood.pdgId == 11, bool),
        #    "is_electron"
        #)
        #xbbVsQCD = self.events.FatJetGood.particleNetMD_Xbb / (self.events.FatJetGood.particleNetMD_Xbb + self.events.FatJetGood.particleNetMD_QCD)
        #self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, xbbVsQCD, 'xbbVsQCD')
        #self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, '2017')

    def process_extra_after_presel(self, variation):
        #FIXME

        #print(self.events.fields)
        
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "btagBB")
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
        if 'DATA' not in self._sample:
            #self.events['passSingleLepElec'] = (ak.count(self.events['ElectronGood']) == 1)
            #self.events['passSingleLepMuon'] = (ak.count(self.events['MuonGood']) == 1)
            ### Add function to implement combinatorics now that we have the sorted list
            self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, self._year)
            print("ele_reco", self.events.ele_reco_sf)
            self.events['ele_id_sf'], self.events['ele_id_sfup'], self.events['ele_id_sfdown'] = sf_ele_id(self.params, self.events, self._year)
            print("ele_id", self.events.ele_id_sf)
            self.events['ele_trig_sf'], self.events['ele_trig_sfup'], self.events['ele_trig_sfdown'] = sf_ele_trig(self.events)
            print("ele_trig", self.events.ele_trig_sf)
            self.events['mu_id_sf'], self.events['mu_id_sfup'], self.events['mu_id_sfdown'] = sf_mu(self.params, self.events, self._year, 'id')
            print("mu_id", self.events.mu_id_sf)
            self.events['mu_iso_sf'], self.events['mu_iso_sfup'], self.events['mu_iso_sfdown'] = sf_mu(self.params, self.events, self._year, 'iso')
            print("mu_iso", self.events.mu_iso_sf)
            self.events['mu_trig_sf'], self.events['mu_trig_sfup'], self.events['mu_trig_sfdown'] = sf_mu(self.params, self.events, self._year, 'trigger')
            print("mu_trig", self.events.mu_trig_sf)
            # NO RECO SF IN RUN 3
            #self.events['mu_reco_sf'], self.events['mu_reco_sfup'], self.events['mu_reco_sfdown'], self.events['mu_reco_sfstat'] = sf_mu(self.params, self.events, self._year, 'reco')
            #print("mu_reco", self.events.mu_reco_sf)
            #corrected_muon_pt, muon_scale_unc = ApplyRochesterCorrections(self._year, self.events.MuonGood, False)
            #self.events['mu_pt_up'], self.events['mu_pt_down'] = corrected_muon_pt + muon_scale_unc, corrected_muon_pt - muon_scale_unc
            #self.events['MuonGood'] = ak.with_field(self.events.MuonGood, corrected_muon_pt, 'pt_corrected')
            #print("Rochester output:", ApplyRochesterCorrections(self._year, self.events.MuonGood, False))
            #btag_variations = ['hf', 'lf', 'hfstats1', 'hfstats2', 'lfstats1', 'lfstats2', 'cferr1', 'cferr2']
            btag_variations = []
            #btag_sf = sf_btag_wp(self.params, self.events.JetGood, self._year, njets=self.events.nJetGood, variations=['central', 'up', 'down'], working_point="M")
            #print(btag_sf['central'][0])
            # btagging
            # 1. Grab the specific lookup tool for this sample chunk
            # If unmapped (or if Data), default to a lambda that returns 1.0s
            eff_lookup = self.btag_eff_maps.get(
                self._sample, 
                lambda pt, eta, flav: ak.ones_like(pt)
            )
            
            # 2. Call the clean wrapper function
            # (Make sure apply_btag_sf is imported at the top of your file!)
            self.events['btag_sf'], self.events['btag_sfup'], self.events['btag_sfdown'] = apply_btag_sf(
                self.params, 
                self.events, 
                self._year, 
                eff_lookup, 
                working_point="M"
            )
            print("BTAG:", self.events['btag_sf'])
            #self.events['btag_sfhf'] = btag_sf['hf'][0]
            #self.events['btag_sfhf_up'] = btag_sf['hf'][1]
            #self.events['btag_sfhf_down'] = btag_sf['hf'][2]
            #self.events['btag_sflf'] = btag_sf['lf'][0]
            #self.events['btag_sflf_up'] = btag_sf['lf'][1]
            #self.events['btag_sflf_down'] = btag_sf['lf'][2]
            # self.events['btag_sfhfstats1'] = btag_sf['hfstats1'][0]
            # self.events['btag_sfhfstats1_up'] = btag_sf['hfstats1'][1]
            # self.events['btag_sfhfstats1_down'] = btag_sf['hfstats1'][2]
            # self.events['btag_sfhfstats2'] = btag_sf['hfstats2'][0]
            # self.events['btag_sfhfstats2_up'] = btag_sf['hfstats2'][1]
            # self.events['btag_sfhfstats2_down'] = btag_sf['hfstats2'][2]
            # self.events['btag_sflfstats1'] = btag_sf['lfstats1'][0]
            # self.events['btag_sflfstats1_up'] = btag_sf['lfstats1'][1]
            # self.events['btag_sflfstats1_down'] = btag_sf['lfstats1'][2]
            # self.events['btag_sflfstats2'] = btag_sf['lfstats2'][0]
            # self.events['btag_sflfstats2_up'] = btag_sf['lfstats2'][1]
            # self.events['btag_sflfstats2_down'] = btag_sf['lfstats2'][2]
            # self.events['btag_sfcferr1'] = btag_sf['cferr1'][0]
            # self.events['btag_sfcferr1_up'] = btag_sf['cferr1'][1]
            # self.events['btag_sfcferr1_down'] = btag_sf['cferr1'][2]
            # self.events['btag_sfcferr2'] = btag_sf['cferr2'][0]
            # self.events['btag_sfcferr2_up'] = btag_sf['cferr2'][1]
            # self.events['btag_sfcferr2_down'] = btag_sf['cferr2'][2]
            #bbtag_sf = sf_bbtag(self.params, self.events.FatJetGood, self._year, njets=self.events.nFatJetGood, variations=['central','up','down'])
            #self.events['bbtag_sf'] = bbtag_sf['central'][0]
            #self.events['bbtag_sfup'] = bbtag_sf['up'][0]
            #self.events['bbtag_sfdown'] = bbtag_sf['down'][0]
            #print("BBTAG SF: ", bbtag_sf.keys())
            # print("BTAG SF: ", len(btag_sf['central'][0]))
            # print("len of events", len(self.events.ele_id_sf))
            # print("number of jets", self.events.nJetGood)
            self.events['puWeight'], self.events['puWeight_up'], self.events['puWeight_down'] = sf_pileup_reweight(self.params, self.events, self._year)
            self.events['isr'], self.events['isr_up'], self.events['isr_down'] = sf_partonshower_isr(self.events)
            self.events['fsr'], self.events['fsr_up'], self.events['fsr_down'] = sf_partonshower_fsr(self.events)
            #self.events['mu_r_up'], self.events['mu_r_down'], self.events['mu_f_up'], self.events['mu_f_down'], self.events['mu_rf_up'], self.events['mu_rf_down'], = sf_scaleweights(self.events)
            #print("fsr_up", self.events.fsr_up)
            #print("isr_up", self.events.isr_up)
            match_gen_lep(self.events)
            #self.do_parton_matching()
            #self.count_partons()
            
        zh_helper(self.events)
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        elif (("TTbb" in self._sample) | ("TTTo" in self._sample) | ("ttbb" in self._sample)):
            match_gen_tt(self.events, self._sample)
        else:
            self.events['process'] = self._sample
            self.events['topptWeight'] = 1
            self.events['topptWeight_Up'] = 1
            self.events['topptWeight_Down'] = 1
        #match_tt_products(self.events)
        #applyDNN(self.events)
        calc_weight(self.events, self.output, self._dataset, self.params, self._year)
        #print("XSEC: ", self.events.metadata['xsec'])
        #print("LUMI: ", self.params.sample_params['lumi']['lumi'])
        #print("genWeights total: ", self.output['sum_signOf_genweights'][self._dataset])
        #if 'TTbb' in self._sample:
        #    add_weights_to_ttbb(self.events, self._sample)

        # Add outputs for data/mc plots
        get_controlvars(self.events)

        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanet1.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetZonly.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalanced24.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedNoTTCC.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedXbbVsQCD.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanet_uscms.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/delete.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/withcuts.onnx")

        self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanet2024_encoderFeatures.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedTruncated.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedAssignment.onnx")

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

    def onnx_inference(self, model_file) -> ak.Array:

        try:
            worker = get_worker()
        except ValueError:
            worker = None
            
        if worker is None:
            import onnxruntime as ort
            sess_options = ort.SessionOptions()
            sess_options.intra_op_num_threads = 1
            sess_options.inter_op_num_threads = 1
            sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            model_session = ort.InferenceSession(
                model_file,
                sess_options = sess_options,
                providers=['CPUExecutionProvider']
            )
        else:
            model_session = worker.data['model_session']

        #print(model_session)

        btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]
        pad_dict = {"btagB":0., "pt":0., "phi":0., "eta":0.}
        jets_padded = ak.zip(
            {key : ak.fill_none(ak.pad_none(self.events.JetGood[key], 16, clip=True), value) for key, value in pad_dict.items()}
        )

        data = np.transpose(
            np.stack([
                np.log(1 + ak.to_numpy(jets_padded.pt)),
                ak.to_numpy(jets_padded.eta),
                np.sin(ak.to_numpy(jets_padded.phi)),
                np.cos(ak.to_numpy(jets_padded.phi)),
                ak.to_numpy(jets_padded.btagB),
            ]),
            axes=[1,2,0]).astype(np.float32)

        mask = ~ak.to_numpy(jets_padded.pt == 0)

        met_data = np.stack([np.log(1+ ak.to_numpy(self.events.PuppiMET.pt)), # pT
                             np.sin(ak.to_numpy(self.events.PuppiMET.phi)), #sin(phi)
                             np.cos(ak.to_numpy(self.events.PuppiMET.phi)) #cos(phi)
                             ], axis=1)[:,None,:].astype(np.float32)

        lep_data = np.stack([np.log(1 + ak.to_numpy(self.events.LeptonGood[:,0].pt)),
                             ak.to_numpy(self.events.LeptonGood[:,0].eta),
                             np.sin(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
                             np.cos(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
                             ], axis=1)[:,None,:].astype(np.float32)

        fatjet_data = np.stack([np.log(1 + ak.to_numpy(self.events.FatJetSorted[:,0].pt)),
                             ak.to_numpy(self.events.FatJetSorted[:,0].eta),
                             np.sin(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
                             np.cos(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
                             ak.to_numpy(self.events.FatJetSorted[:,0].btagBB),
                             ], axis=1)[:,None,:].astype(np.float32)
        event_data = ak.to_numpy(self.events.n_b_inZH[:,None,None]).astype(np.float32)

        mask_global = np.ones(shape=[met_data.shape[0], 1]) == 1
        if "Zonly" in model_file:
            output_names = ["EVENT/ttzbb"]
        elif "1" in model_file:
            output_names = ["EVENT/tthbb", "EVENT/ttbb", "EVENT/ttcc", "EVENT/ttlf"]
        #else:
        #output_names = ["EVENT/ttzbb", "EVENT/tthbb", "EVENT/ttbb", "EVENT/ttcc", "EVENT/ttlf"]
        #output_names = ["EVENT/ttzbb", "EVENT/tthbb", "EVENT/ttbb", "EVENT/ttlf"]
        output_names = ["EVENT/signal", "encoder_event_vector"]
        outputs = model_session.run(input_feed={
            "Jet_data": data,
            "Jet_mask": mask,
            "Met_data": met_data,
            "Met_mask": mask_global,
            "Lepton_data": lep_data,
            "Lepton_mask": mask_global,
            "FatJet_data": fatjet_data,
            "FatJet_mask": mask_global},
            #"Event_data": event_data,
            #"Event_mask": mask_global},
        output_names=output_names
        )
        #print("-" * 60)
        #print("OUTPUTS")
        #print(outputs)
        #print("-" * 60)
        #print((ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])))
       

        outputs_zipped = dict(zip(output_names, outputs))
        print(outputs_zipped.keys())
        keys = list(outputs_zipped.keys())
        #print((ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) for key, value in outputs_zipped.items())
        #print("0", (ak.from_numpy(value[:,0]) for key, value in outputs_zipped.items()))
        #print("1", (ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()))
        #print("2", (ak.from_numpy(value[:,2]) for key, value in outputs_zipped.items()))
        #print("3", (ak.from_numpy(value[:,3]) for key, value in outputs_zipped.items()))
        #print("4", (ak.from_numpy(value[:,4]) for key, value in outputs_zipped.items()))
        print(len(outputs_zipped[keys[1]][0,:]))
        if "Zonly" in model_file:
            self.events["spanet_outputZ"] = ak.zip(
                {
                    key.split("/")[-1]: ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()
                }
            )
        #else:
        self.events["spanet_output"] = ak.zip(
            {
                #key.split("/")[-1]:  (ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) for key, value in outputs_zipped.items()
                #keys[0].split("/")[-1]:  (ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2])) / (ak.from_numpy(outputs_zipped[keys[0]][:,0]) + ak.from_numpy(outputs_zipped[keys[0]][:,3]) + ak.from_numpy(outputs_zipped[keys[0]][:,4]) + ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2]))
                keys[0].split("/")[-1]:  ak.from_numpy(outputs_zipped[keys[0]][:,1]) / (ak.from_numpy(outputs_zipped[keys[0]][:,0]) + ak.from_numpy(outputs_zipped[keys[0]][:,1]) + ak.from_numpy(outputs_zipped[keys[0]][:,2]) + ak.from_numpy(outputs_zipped[keys[0]][:,3]))
                #"out_1": ak.from_numpy(value[:,1] for key, value in outputs_zipped.items()),
                #"out_2": ak.from_numpy(value[:,2] for key, value in outputs_zipped.items()),
            }
        )
        self.events["sig_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,1])
        self.events["ttbb_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,2])
        self.events["ttlf_score"] = ak.from_numpy(outputs_zipped[keys[0]][:,3])
        self.events["encoderFeature1"] = ak.from_numpy(outputs_zipped[keys[1]][:,0])
        self.events["encoderFeature2"] = ak.from_numpy(outputs_zipped[keys[1]][:,1])
        self.events["encoderFeature3"] = ak.from_numpy(outputs_zipped[keys[1]][:,2])
        self.events["encoderFeature4"] = ak.from_numpy(outputs_zipped[keys[1]][:,3])
        self.events["encoderFeature5"] = ak.from_numpy(outputs_zipped[keys[1]][:,4])
        self.events["encoderFeature6"] = ak.from_numpy(outputs_zipped[keys[1]][:,5])
        self.events["encoderFeature7"] = ak.from_numpy(outputs_zipped[keys[1]][:,6])
        self.events["encoderFeature8"] = ak.from_numpy(outputs_zipped[keys[1]][:,7])
        self.events["encoderFeature9"] = ak.from_numpy(outputs_zipped[keys[1]][:,8])
        self.events["encoderFeature10"] = ak.from_numpy(outputs_zipped[keys[1]][:,9])
        self.events["encoderFeature11"] = ak.from_numpy(outputs_zipped[keys[1]][:,10])
        self.events["encoderFeature12"] = ak.from_numpy(outputs_zipped[keys[1]][:,11])
        self.events["encoderFeature13"] = ak.from_numpy(outputs_zipped[keys[1]][:,12])
        self.events["encoderFeature14"] = ak.from_numpy(outputs_zipped[keys[1]][:,13])
        self.events["encoderFeature15"] = ak.from_numpy(outputs_zipped[keys[1]][:,14])
        self.events["encoderFeature16"] = ak.from_numpy(outputs_zipped[keys[1]][:,15])
        self.events["encoderFeature17"] = ak.from_numpy(outputs_zipped[keys[1]][:,16])
        self.events["encoderFeature18"] = ak.from_numpy(outputs_zipped[keys[1]][:,17])
        self.events["encoderFeature19"] = ak.from_numpy(outputs_zipped[keys[1]][:,18])
        self.events["encoderFeature20"] = ak.from_numpy(outputs_zipped[keys[1]][:,19])
        self.events["encoderFeature21"] = ak.from_numpy(outputs_zipped[keys[1]][:,20])
        self.events["encoderFeature22"] = ak.from_numpy(outputs_zipped[keys[1]][:,21])
        self.events["encoderFeature23"] = ak.from_numpy(outputs_zipped[keys[1]][:,22])
        self.events["encoderFeature24"] = ak.from_numpy(outputs_zipped[keys[1]][:,23])
        self.events["encoderFeature25"] = ak.from_numpy(outputs_zipped[keys[1]][:,24])
        self.events["encoderFeature26"] = ak.from_numpy(outputs_zipped[keys[1]][:,25])
        self.events["encoderFeature27"] = ak.from_numpy(outputs_zipped[keys[1]][:,26])
        self.events["encoderFeature28"] = ak.from_numpy(outputs_zipped[keys[1]][:,27])
        self.events["encoderFeature29"] = ak.from_numpy(outputs_zipped[keys[1]][:,28])
        self.events["encoderFeature30"] = ak.from_numpy(outputs_zipped[keys[1]][:,29])
        self.events["encoderFeature31"] = ak.from_numpy(outputs_zipped[keys[1]][:,30])
        self.events["encoderFeature32"] = ak.from_numpy(outputs_zipped[keys[1]][:,31])
        print(self.events.encoderFeature20)
        #print(self.events.spanet_outputH)
