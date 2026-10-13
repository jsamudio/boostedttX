import awkward as ak
import numpy as np
import correctionlib

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
#from pocket_coffea.lib.weights_manager import WeightsManager
from object_cleaning_functions import soft_lep_sel, lep_sel, bjet_sel, qjet_sel, lep_softlep_combo, jet_selection_custom
from custom_cut_functions import sortbyscore
from analysis_years import xbb_wp_cuts
from cand_helper import zh_helper#, ak4_truncate
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig
#import dnn_model
#from applyDNN import applyDNN
from weight_handler import calc_weight, add_weights_to_ttbb
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from custom_scale_factors import sf_mu, sf_partonshower_isr, sf_partonshower_fsr, sf_scaleweights, apply_jms_jmr

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

class ZHbbBaseProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        self.pnet_mass_cset = correctionlib.CorrectionSet.from_file(self.workflow_options["pnet_mass_json"])

    def process_extra_before_skim(self):
        self.events['sum_sign_genw'] = np.sum(np.sign(self.events['genWeight']))


    def apply_object_preselection(self, variation):
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
        print("quarks", ak.num(quarks, axis=-1))

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
            print("Matched fatjets:", matched_fatjets)

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
        print("quarks", quarks.pdgId)
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
        print("Higgs:", higgs.provenance)
        #self.events["HiggsMatched"] = higgs

        # Calling our general object_matching function.
        # The output is an awkward array with the shape of the second argument and None where there is no matching.
        # So, calling like this, we will get out an array of matched_quarks with the dimension of the JetGood.
        # FIXME hardcoded truncation of AK4 jets for the matched jet collection
        matched_quarks, matched_jets, deltaR_matched = object_matching(
            #quarks, self.events.JetGoodTruncated, dr_min=self.dr_min
            quarks, self.events.JetGood, dr_min=self.dr_min
        )     
        matched_higgs, matched_fatjets, deltaR_matchedAK8 = object_matching(
            higgs, self.events.FatJetSorted, dr_min=0.6
        )
        self.events["FatJetMatched"] = ak.with_field(
            matched_fatjets, deltaR_matchedAK8, "dRMatchedJet"
        )
        print("Higgs:", matched_higgs)
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

    def define_common_variables_after_presel(self, variation):
        super().define_common_variables_before_presel(variation=variation)

        # Define labels for btagged jets at different working points
        for wp, val in self.params.btagging.working_point[self._year]["btagging_WP"].items():
            self.events["JetGood"] = ak.with_field(
                self.events.JetGood,
                ak.values_astype(self.events.JetGood.btagDeepFlavB > val, int),
                f"btag_{wp}"
            )
        #self.events.FatJetGood['xbbVsQCD'] = self.events.FatJetGood.particleNetMD_Xbb / (self.events.FatJetGood.particleNetMD_Xbb + self.events.FatJetGood.particleNetMD_QCD)
        xbbVsQCD = self.events.FatJetGood.particleNetMD_Xbb / (self.events.FatJetGood.particleNetMD_Xbb + self.events.FatJetGood.particleNetMD_QCD)
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, xbbVsQCD, 'xbbVsQCD')

        # Standard ParticleNet-MD XbbvsQCD Run 2 UL cuts (analysis_years.XBB_WP_DICT)
        wp_cuts = xbb_wp_cuts(self._year)

        # Add boolean flags for the particleNetMD_XbbvsQCD Working Points
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, ak.values_astype(xbbVsQCD > wp_cuts["L"], int), 'xbb_wp_L')
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, ak.values_astype(xbbVsQCD > wp_cuts["M"], int), 'xbb_wp_M')
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, ak.values_astype(xbbVsQCD > wp_cuts["T"], int), 'xbb_wp_T')


    def process_extra_after_presel(self, variation):
        print(self.events.FatJetGood.xbbVsQCD)
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")

        # --- AK4 / AK8 Overlap Removal ---
        # Get the Z/H candidate (the leading FatJet sorted by Xbb score)
        leading_fatjet = ak.firsts(self.events.FatJetSorted)
        # Compute deltaR between all AK4 JetGood and the leading AK8
        dr_jets = self.events.JetGood.delta_r(leading_fatjet)
        # Keep AK4 jets with dR > 0.8 (standard AK8 cone size) or events with no AK8
        overlap_mask = (dr_jets > 0.8) | ak.is_none(dr_jets)
        self.events["JetGood"] = self.events.JetGood[overlap_mask]
        
        # Re-apply bJet and qJet selections just in case they are used downstream by zh_helper
        self.events['bJetGood'] = bjet_sel(self.events.JetGood, self.params, self._year)
        self.events['qJetGood'] = qjet_sel(self.events.JetGood, self.params, self._year)
        # ---------------------------------
        
        #self.events['passSingleLepElec'] = (ak.count(self.events['ElectronGood']) == 1)
        #self.events['passSingleLepMuon'] = (ak.count(self.events['MuonGood']) == 1)
        ### Add function to implement combinatorics now that we have the sorted list
        zh_helper(self.events)
        #ak4_truncate(self.events)
        match_gen_lep(self.events)
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        else:
            match_gen_tt(self.events, self._sample)
        #applyDNN(self.events)
        calc_weight(self.events, self.output, self._dataset, self.params, self._year)
        print("XSEC: ", self.events.metadata['xsec'])
        print("LUMI: ", self.params.lumi['picobarns'][self._year]['tot'])
        print("genWeights total: ", self.output['sum_signOf_genweights'][self._dataset])
        #if 'TTbb' in self._sample:
        #    add_weights_to_ttbb(self.events, self._sample)
        self.do_parton_matching()
        self.count_partons()

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



