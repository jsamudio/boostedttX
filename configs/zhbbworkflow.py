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

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

class ZHbbBaseProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

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



