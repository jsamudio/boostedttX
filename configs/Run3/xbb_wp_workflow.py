"""Workflow of the Xbb working-point derivation (Run 3): all fatjets of the analysis preselection,
with their flavour (B/C hadrons) and the match to a hard-process H/Z -> bb, for derive_xbb_wps.py."""
import awkward as ak

from pocket_coffea.utils.configurator import Configurator
from custom_cut_functions import sortbyscore
from genmatcher import match_gen_lep
from parton_matching import hard_bosons, hard_partons, match_fatjets
from spanet_inference import ZHbbRun3Processor


class ZHbbRun3XbbWPProcessor(ZHbbRun3Processor):
    '''
    Same objects and preselection as the analysis (no b-tag category, no Xbb requirement).
    Adds to FatJetGood (MC):
      matched      - a hard-process H/Z (-> quarks) within dR < 0.6
      n_boson_b    - b quarks from the H/Z within dR < 0.8 (2 = fully contained H/Z -> bb)
    and the process labels / gen flags used by the subsamples (tt+LF, Zbb).
    '''

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "btagBB")
        if not self._isMC:
            return
        match_gen_lep(self.events)
        self.label_process()

        fatjets = match_fatjets(self.events.FatJetGood, hard_bosons(self.events.GenPart), 0.6)
        partons = hard_partons(self.events.GenPart)
        boson_b = partons[(partons.provenance == 1) & (abs(partons.pdgId) == 5)]
        dr = fatjets.metric_table(boson_b)
        fatjets = ak.with_field(fatjets, ak.sum(dr < 0.8, axis=2), "n_boson_b")
        self.events['FatJetGood'] = fatjets
