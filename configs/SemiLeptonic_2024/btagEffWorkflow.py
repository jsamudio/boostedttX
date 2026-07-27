import awkward as ak
import numpy as np
import hist
from pocket_coffea.workflows.base import BaseProcessorABC
from object_cleaning_functions import soft_lep_sel, lep_sel, bjet_sel, qjet_sel, lep_softlep_combo, jet_selection_custom
from custom_cut_functions import sortbyscore

class BTagEffProcessor(BaseProcessorABC):
    def __init__(self, cfg):
        super().__init__(cfg)
        
        # Define the custom histogram format for the accumulator
        # We use pT and absolute eta to map the efficiencies
        self.output_format["btag_eff"] = hist.Hist(
            hist.axis.StrCategory(["all", "pass"], name="passWP"),
            hist.axis.IntCategory([0, 4, 5], name="flavor", growth=True),
            hist.axis.Variable([20, 30, 50, 70, 100, 140, 200, 300, 600, 1000], name="pt"),
            hist.axis.Regular(5, 0.0, 2.5, name="abseta")
        )

    def apply_object_preselection(self, variation):
        # EXACT copy of your preselection logic so the phase space matches perfectly
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(self.events.Electron, electron_etaSC, "etaSC")

        self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        self.events['SoftMuonGood'] = soft_lep_sel(self.events, "Muon", self.params)
        self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)
        self.events['SoftElectronGood'] = soft_lep_sel(self.events, "Electron", self.params)

        leptons = ak.with_name(
                ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis = 1),
                name='PtEtaPhiMCandidate')
        self.events['LeptonGood'] = leptons[ak.argsort(leptons.pt, ascending=False)]

        self.events['JetGood'], self.jetGoodMask = jet_selection_custom(
            self.events, "Jet", self.params, self._year, "LeptonGood", 
            self.params.object_preselection.btag.tagger
        ) 
        self.events['FatJetGood'], self.fatJetGoodMask = jet_selection_custom(self.events, "FatJet", self.params, self._year, "LeptonGood", "GloParT") # they use some MSD cut which we don't want
        self.events['e_softe'] = lep_softlep_combo(self.events.ElectronGood, self.events.SoftElectronGood)
        self.events['mu_softmu'] = lep_softlep_combo(self.events.MuonGood, self.events.SoftMuonGood)

        # 2. Add the custom fields for the Efficiency axes
        wp_val = self.params.btagging.working_point[self._year]["btagging_WP"]["M"]
        btag_algo = self.params.btagging.working_point[self._year]["btagging_algorithm"]
        
        # Create a boolean mask for jets that passed
        pass_mask = self.events.JetGood[btag_algo] >= wp_val
        
        # Create the sub-collection of passing jets
        self.events["JetGood_Pass"] = self.events.JetGood[pass_mask]
        
        # Clean the flavor to strictly 0, 4, 5 and calculate absolute eta
        for coll in ["JetGood", "JetGood_Pass"]:
            raw_flav = self.events[coll].hadronFlavour
            clean_flav = ak.where((raw_flav != 5) & (raw_flav != 4), 0, raw_flav)
            
            self.events[coll] = ak.with_field(self.events[coll], clean_flav, "clean_flavor")
            self.events[coll] = ak.with_field(self.events[coll], np.abs(self.events[coll].eta), "abseta")
        
    def count_objects(self, variation):
        self.events['nMuonGood'] = ak.num(self.events.MuonGood)
        self.events['nElectronGood'] = ak.num(self.events.ElectronGood)
        self.events['nSoftMuonGood'] = ak.num(self.events.SoftMuonGood)
        self.events['nSoftElectronGood'] = ak.num(self.events.SoftElectronGood)
        self.events['nJetGood'] = ak.num(self.events.JetGood)
        self.events['nFatJetGood'] = ak.num(self.events.FatJetGood)
        self.events['nLeptonGood'] = (self.events['nMuonGood'] + self.events['nElectronGood'])
        
    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "btagBB")
        ZHCand = self.events.FatJetSorted[:,0]
        self.events['ZH_pt'] = ZHCand.pt
        self.events['ZH_M'] = (ZHCand.globalParT3_massCorrX2p * ZHCand.mass * (1-ZHCand.rawFactor))
        self.events['ZH_xbb'] = ZHCand.btagBB        