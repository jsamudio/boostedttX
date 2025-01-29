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
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig
import dnn_model
from applyDNN import applyDNN
from weight_handler import calc_weight, add_weights_to_ttbb
from coffea.analysis_tools import PackedSelection
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

class ZHbbBaseProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)

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

    def do_parton_matching_ttHbb(self) -> ak.Array:

        genparts = self.events.GenPart
        
        children_idxG = ak.without_parameters(genparts.childrenIdxG, behavior={})
        children_idxG_flat = ak.flatten(children_idxG, axis=1)
        genparts_flat = ak.flatten(genparts)
        genparts_offsets = np.concatenate([[0],np.cumsum(ak.to_numpy(ak.num(genparts, axis=1), allow_missing=True))])
        local_index_all = ak.local_index(genparts, axis=1)

        # Get the initial partons, first copy
        initial = genparts.genPartIdxMother == 0
        hard_process = (genparts.hasFlags(['fromHardProcess','isPrompt','isHardProcess', 'isFirstCopy'])) & (genparts.status !=21) # exclude incoming particles
        genparts_hard = genparts[hard_process]

        
        higgs = genparts_hard[genparts_hard.pdgId == 25]
        top = genparts_hard[genparts_hard.pdgId == 6]
        antitop = genparts_hard[genparts_hard.pdgId == -6]
        
        # I want to take hardProcess, final state, BEFORE FSR particles, which have higgs, top, antitop as parents
        # These will become the hard particles 
        from_higgs = genparts_hard.parent.pdgId == 25
        from_top = genparts_hard.parent.pdgId == 6
        from_antitop = genparts_hard.parent.pdgId == -6

        genparts_initial = genparts[initial]

        isr = genparts_initial[(genparts_initial.hasFlags(['fromHardProcess','isPrompt','isHardProcess'])) &
                                  (genparts_initial.status == 23) & (genparts_initial.pdgId != 25)& (abs(genparts_initial.pdgId) != 6)]
        has_isr = ak.num(isr)!=0
        isr_idx = ak.flatten(ak.fill_none(ak.pad_none(local_index_all[initial][(genparts_initial.hasFlags(['fromHardProcess','isPrompt','isHardProcess'])) &
                                  (genparts_initial.status == 23) & (genparts_initial.pdgId != 25)& (abs(genparts_initial.pdgId) != 6)], 1), 0))

        ######
        # Get the hard process particles
        part_from_top = genparts_hard[from_top]
        W_from_top = ak.flatten(part_from_top[part_from_top.pdgId ==24])
        b_from_top = ak.flatten(part_from_top[part_from_top.pdgId ==5])
        part_from_antitop = genparts_hard[from_antitop]
        W_from_antitop = ak.flatten(part_from_antitop[part_from_antitop.pdgId == -24])
        b_from_antitop = ak.flatten(part_from_antitop[part_from_antitop.pdgId == -5])

        b_from_top_idx = ak.flatten(local_index_all[hard_process][from_top][part_from_top.pdgId ==5])
        W_from_top_idx = ak.flatten(local_index_all[hard_process][from_top][part_from_top.pdgId == 24])
        b_from_antitop_idx = ak.flatten(local_index_all[hard_process][from_antitop][part_from_antitop.pdgId == -5])
        W_from_antitop_idx = ak.flatten(local_index_all[hard_process][from_antitop][part_from_antitop.pdgId == -24])

        # This works because they are the firstCopy of the hard_process particles with higgs as parent. We are skipping all the decay chain of the higgs
        b_from_higgs_idx = local_index_all[hard_process][from_higgs]
        b_from_higgs = genparts_hard[from_higgs]

        #--------------
        # Converting local index to global index, already corrected with the offsets
        b_from_top_idxG = ak.to_numpy(b_from_top_idx + genparts_offsets[:-1], allow_missing=False)
        W_from_top_idxG = ak.to_numpy(W_from_top_idx + genparts_offsets[:-1], allow_missing=False)
        b_from_antitop_idxG = ak.to_numpy(b_from_antitop_idx + genparts_offsets[:-1], allow_missing=False)
        W_from_antitop_idxG = ak.to_numpy(W_from_antitop_idx + genparts_offsets[:-1], allow_missing=False)
        isr_idxG = ak.to_numpy(isr_idx + genparts_offsets[:-1], allow_missing=False)
        b_from_higgs_idxG = ak.to_numpy(b_from_higgs_idx + genparts_offsets[:-1], allow_missing=False)

        # Some inputs needed for numba functions for indexing
        
        genparts_flat_eta = ak.without_parameters(genparts_flat.eta, behavior={})
        genparts_flat_phi = ak.without_parameters(genparts_flat.phi, behavior={})
        genparts_flat_pt = ak.without_parameters(genparts_flat.pt, behavior={})
        genparts_flat_pdgId = ak.without_parameters(genparts_flat.pdgId, behavior={})
        genparts_flat_statusFlags = ak.without_parameters(genparts_flat.statusFlags, behavior={})
        
        firstgenpart_idxG = ak.firsts(genparts[:,0].children).genPartIdxMotherG
        firstgenpart_idxG_numpy = ak.to_numpy( firstgenpart_idxG, allow_missing=False)
        local_ind = ak.to_numpy(ak.local_index(firstgenpart_idxG), allow_missing=False)
        nevents = firstgenpart_idxG_numpy.shape[0]

        ### Analyze the W decay
        W_from_top_islep, W_from_top_decay = analyze_W_flat( W_from_top_idxG, 
                                                             children_idxG_flat,
                                                             genparts_flat_statusFlags,
                                                             genparts_flat_pdgId,
                                                             firstgenpart_idxG_numpy,
                                                             genparts_offsets,
                                                             nevents)
                                                            

        W_from_antitop_islep, W_from_antitop_decay = analyze_W_flat( W_from_antitop_idxG, 
                                                                     children_idxG_flat,
                                                                     genparts_flat_statusFlags,
                                                                     genparts_flat_pdgId,
                                                                     firstgenpart_idxG_numpy,
                                                                     genparts_offsets,
                                                                     nevents)
        
        # assuming semilep only
        W_had_decay_idx = np.where(~W_from_top_islep[:,None],W_from_top_decay, W_from_antitop_decay )
        W_lep_decay_idx = np.where(W_from_top_islep[:,None], W_from_top_decay, W_from_antitop_decay )

        # Now getting all the global Idx of particles for which we want to analyze the decays chain,
        # looking for the highest pt emission after FSR radiation
        part_input_G  = np.concatenate([b_from_top_idxG[:,None],
                                        b_from_antitop_idxG[:,None],
                                        b_from_higgs_idxG,
                                        isr_idxG[:,None],
                                        W_had_decay_idx,
                                        ], axis=1)

        parton_decay_id = analyze_parton_decays_flat_nomesons(part_input_G,
                                                              children_idxG_flat,
                                                              genparts_flat_eta,
                                                              genparts_flat_phi,
                                                              genparts_flat_pt,
                                                              genparts_flat_pdgId,
                                                              self.dr_min_postfsr,
                                                              firstgenpart_idxG_numpy,
                                                              genparts_offsets,
                                                              nevents)

        
        b_from_top_lastcopy = genparts_flat[parton_decay_id[:,0]]
        b_from_antitop_lastcopy = genparts_flat[parton_decay_id[:,1]]
        b_from_higgs_lastcopy = genparts_flat[parton_decay_id[:,2:4]]
        isr_lastcopy = genparts_flat[parton_decay_id[:,4]]
        part_from_Whad_lastcopy = genparts_flat[parton_decay_id[:,5]]
        # We don't take the last copy of the leptonic particles, but the born one
        part_from_Wlep = genparts_flat[W_lep_decay_idx]
        
        # Now we can perform the genmatching
        quarks_initial = genparts_flat[part_input_G]
        quarks_lastcopy = genparts_flat[parton_decay_id]

        quarks_provenance = np.zeros(parton_decay_id.shape, dtype=np.int32)
        quarks_provenance[:, 0] = np.where(W_from_top_islep, 3, 2)  #b from leptonic top =3, hadronic =2
        quarks_provenance[:, 1] = np.where(W_from_antitop_islep, 3, 2)
        quarks_provenance[:, 2:4] = 1 #Higgs
        quarks_provenance[:, 4] = 4 #isr
        quarks_provenance[:, 5:7] = 5 #W hadronic decay

        # Assign provenance
        # 1 - from higgs
        # 2 - from top hadronic bquark
        # 3 - from top leptonic bquark
        # 4 - from ISR
        # 5 - from W hadronic decay
        quarks_initial["provenance"] = ak.Array(quarks_provenance, behavior={})
        quarks_lastcopy["provenance"] = quarks_initial["provenance"]
        

        # Calling our general object_matching function.
        # The output is an awkward array with the shape of the second argument and None where there is no matching.
        # So, calling like this, we will get out an array of matched_quarks with the dimension of the JetGood.
        matched_quarks, matched_jets, deltaR_matched = object_matching(
            quarks_lastcopy, self.events.JetGood, dr_min=self.dr_min
        )

        #Saving stuff
        self.events["JetGoodMatched"] = ak.with_field(
            matched_jets, deltaR_matched, "dRMatchedJet"
        )
        self.events["JetGoodMatched"] = ak.with_field(
            self.events.JetGoodMatched, matched_quarks.provenance, "provenance")
        
        self.events["PartonInitial"] = quarks_initial
        self.events["PartonLastCopy"] = quarks_lastcopy
        # Saving the matched partons only
        self.events["PartonLastCopyMatched"] = matched_quarks
        self.matched_partons_mask = ~ak.is_none(self.events.JetGoodMatched, axis=1)

        self.events["LeptonGenLevel"] = part_from_Wlep
        self.events["HiggsGen"] = higgs
        
        self.events["TopGen"] = top
        self.events["AntiTopGen"] = antitop
        self.events["TopGen_islep"] = W_from_top_islep
        self.events["AntiTopGen_islep"] = W_from_antitop_islep
        

    def count_partons(self):
        self.events["nPartonLastCopyMatched"] = ak.count(
            self.events.PartonLastCopyMatched.pt, axis=1
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


    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "particleNetMD_Xbb")
        #self.events['passSingleLepElec'] = (ak.count(self.events['ElectronGood']) == 1)
        #self.events['passSingleLepMuon'] = (ak.count(self.events['MuonGood']) == 1)
        ### Add function to implement combinatorics now that we have the sorted list
        zh_helper(self.events)
        match_gen_lep(self.events)
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        else:
            match_gen_tt(self.events, self._sample)
        applyDNN(self.events)
        calc_weight(self.events, self.output, self._dataset, self.params)
        print("XSEC: ", self.events.metadata['xsec'])
        print("LUMI: ", self.params.sample_params['lumi']['lumi'])
        print("genWeights total: ", self.output['sum_signOf_genweights'][self._dataset])
        if 'TTbb' in self._sample:
            add_weights_to_ttbb(self.events, self._sample)
        if self._sample in ["ttHTobb"]:
            self.do_parton_matching_ttHbb()
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



