import awkward as ak
import numpy as np

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
#from pocket_coffea.lib.weights_manager import WeightsManager
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    jet_selection,
    btagging,
    get_dilepton,
)
from object_cleaning_functions import soft_lep_sel, lep_sel, fatjet_sel, bjet_sel, qjet_sel, lep_softlep_combo, jet_sel, fatjet_sel2
from custom_cut_functions import sortbyscore
from cand_helper import zh_helper
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig, match_tt_products
import dnn_model
from applyDNN import applyDNN
from weight_handler import calc_weight, add_weights_to_ttbb
from coffea.analysis_tools import PackedSelection
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.scale_factors import sf_ele_reco, sf_ele_id, sf_mu, sf_btag
from dask.distributed import get_worker

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

class ZHbbSpanetProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        if not "spanet_model" in self.workflow_options:
            raise ValueError("Key `spanet_model` not found in workflow options. Please specify the path to the ONNX model.")
        elif not self.workflow_options["spanet_model"].endswith(".onnx"):
            raise ValueError("Key `spanet_model` should be the path of an ONNX model.")
        

    def process_extra_before_skim(self):
        self.events['sum_sign_genw'] = np.sum(np.sign(self.events['genWeight']))

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

        self.events['JetGood'], self.jetGoodMask = jet_sel(self.events, "Jet", self.params, "LeptonGood")
        self.events['FatJetGood'] = fatjet_sel(self.events, self.params, "LeptonGood")
        self.events['FatJetGood2'] = fatjet_sel2(self.events, self.params, "LeptonGood")
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
        matched_quarks, matched_jets, deltaR_matched = object_matching(
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
        #self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, '2017')


    def process_extra_after_presel(self, variation):
        #FIXME
        #self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "particleNetMD_Xbb")
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")
        #self.events['passSingleLepElec'] = (ak.count(self.events['ElectronGood']) == 1)
        #self.events['passSingleLepMuon'] = (ak.count(self.events['MuonGood']) == 1)
        ### Add function to implement combinatorics now that we have the sorted list
        self.events['ele_reco_sf'], self.events['ele_reco_sfup'], self.events['ele_reco_sfdown'] = sf_ele_reco(self.params, self.events, self._year)
        self.events['ele_id_sf'], self.events['ele_id_sfup'], self.events['ele_id_sfdown'] = sf_ele_reco(self.params, self.events, self._year)
        self.events['mu_id_sf'], self.events['mu_id_sfup'], self.events['mu_id_sfdown'] = sf_mu(self.params, self.events, self._year, 'id')
        self.events['mu_iso_sf'], self.events['mu_iso_sfup'], self.events['mu_iso_sfdown'] = sf_mu(self.params, self.events, self._year, 'iso')
        btag_variations = self.params.systematic_variations.weight_variations.sf_btag[self._year]
        btag_sf = sf_btag(self.params, self.events.JetGood, self._year, njets=self.events.nJetGood, variations=['central']+btag_variations)
        #bbtag_sf = sf_btag(self.params, self.event.FatJetGood, self._year, njets=self.events.nFatJetGood, 
        print("BTAG SF: ", btag_sf.keys())
        zh_helper(self.events)
        match_gen_lep(self.events)
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        else:
            match_gen_tt(self.events, self._sample)
            #self.get_ttbb_LHE_info()
        #match_tt_products(self.events)
        #applyDNN(self.events)
        calc_weight(self.events, self.output, self._dataset, self.params)
        print("XSEC: ", self.events.metadata['xsec'])
        print("LUMI: ", self.params.sample_params['lumi']['lumi'])
        print("genWeights total: ", self.output['sum_signOf_genweights'][self._dataset])
        if 'TTbb' in self._sample:
            add_weights_to_ttbb(self.events, self._sample)
        self.do_parton_matching()
        self.count_partons()
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanet1.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetZonly.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalanced24.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedNoTTCC.onnx")
        #self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalancedXbbVsQCD.onnx")
        self.onnx_inference(model_file=f"/cms/data/jsamudio/boosted/boostedttX/configs/spanet_uscms.onnx")
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
            sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            model_session = ort.InferenceSession(
                model_file,
                sess_options = sess_options,
                providers=['CPUExecutionProvider']
            )
        else:
            model_session = worker.data['model_session']

        print(model_session)

        btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]
        pad_dict = {btagging_algorithm:0., "pt":0., "phi":0., "eta":0.}
        jets_padded = ak.zip(
            {key : ak.fill_none(ak.pad_none(self.events.JetGood[key], 16, clip=True), value) for key, value in pad_dict.items()}
        )

        data = np.transpose(
            np.stack([
                np.log(1 + ak.to_numpy(jets_padded.pt)),
                ak.to_numpy(jets_padded.eta),
                np.sin(ak.to_numpy(jets_padded.phi)),
                np.cos(ak.to_numpy(jets_padded.phi)),
                ak.to_numpy(jets_padded.btagDeepFlavB),
            ]),
            axes=[1,2,0]).astype(np.float32)

        mask = ~ak.to_numpy(jets_padded.pt == 0)

        met_data = np.stack([np.log(1+ ak.to_numpy(self.events.MET.pt)),
                             ak.zeros_like(self.events.MET.pt).to_numpy(),
                             np.sin(ak.to_numpy(self.events.MET.phi)),
                             np.cos(ak.to_numpy(self.events.MET.phi))
                             ], axis=1)[:,None,:].astype(np.float32)

        lep_data = np.stack([np.log(1 + ak.to_numpy(self.events.LeptonGood[:,0].pt)),
                             ak.to_numpy(self.events.LeptonGood[:,0].eta),
                             np.sin(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
                             np.cos(ak.to_numpy(self.events.LeptonGood[:,0].phi)),
                             ak.to_numpy(self.events.LeptonGood[:,0].is_electron).astype(np.int32),
                             ], axis=1)[:,None,:].astype(np.float32)

        fatjet_data = np.stack([np.log(1 + ak.to_numpy(self.events.FatJetSorted[:,0].pt)),
                             ak.to_numpy(self.events.FatJetSorted[:,0].eta),
                             np.sin(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
                             np.cos(ak.to_numpy(self.events.FatJetSorted[:,0].phi)),
                             ak.to_numpy(self.events.FatJetSorted[:,0].xbbVsQCD),
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
        output_names = ["EVENT/signal"]
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
        #print(outputs)
        #print(ak.from_numpy(outputs[:,1]) + ak.from_numpy(outputs[:,2])) / (ak.from_numpy(outputs[:,3]) + ak.from_numpy(outputs[:,4]) + ak.from_numpy(outputs[:,1]) + ak.from_numpy(outputs[:,2]))

        outputs_zipped = dict(zip(output_names, outputs))
        #print(outputs_zipped)
        if "Zonly" in model_file:
            self.events["spanet_outputZ"] = ak.zip(
                {
                    key.split("/")[-1]: ak.from_numpy(value[:,1]) for key, value in outputs_zipped.items()
                }
            )
        #else:
        self.events["spanet_outputH"] = ak.zip(
            {
                key.split("/")[-1]:  (ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) / (ak.from_numpy(value[:,0]) + ak.from_numpy(value[:,3]) + ak.from_numpy(value[:,4]) + ak.from_numpy(value[:,1]) + ak.from_numpy(value[:,2])) for key, value in outputs_zipped.items()
            }
        )
        print(self.events.spanet_outputH)
