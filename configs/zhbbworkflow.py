import awkward as ak

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    jet_selection,
    btagging,
    get_dilepton,
)
from custom_cut_functions import calc_cutBasedNoIso, lep_sel

class ZHbbBaseProcessor (BaseProcessorABC):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

## Perform precut during skim? I think skim is specififed in config outside of base processor
    '''
    @staticmethod
    def is_lep_cleaned(lep, l_k, jets, j_k, cut=0.4): # cut should be 0.4, 0.8
        # cleans any jets matched with analysis lepton
        lep_eta, lep_phi = lep[f'{l_k}_eta'], lep[f'{l_k}_phi']
        jets_eta, jets_phi = jets[f'{j_k}_eta'], jets[f'{j_k}_phi']
        jet_mask = ak_crosscleaned(lep_eta,lep_phi,jets_eta,jets_phi,cut)
        return jet_mask
    '''

    #def define_common_variables_before_presel(self, variation):
        #self.events["Electron_cutBasedNoIso"] = calc_cutBasedNoIso(self.events, self.params)

    def apply_object_preselection(self, variation):
        #soft e, e, soft mu, mu, jet, fatjet
        #self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        self.events['MuonGood'] = self.events.Muon

        #self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)

    def count_objects(self, variation):
        self.events['nMuonGood'] = ak.num(self.events.MuonGood)


