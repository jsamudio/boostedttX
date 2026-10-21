"""Workflow of the SPANet training inputs (Run 2): the analysis objects + gen matching, no SFs."""
import awkward as ak
import vector

from pocket_coffea.utils.configurator import Configurator
from custom_cut_functions import sortbyscore, get_spanet_selection
from cand_helper import zh_helper
from genmatcher import match_gen_lep
from parton_matching import do_parton_matching
from event_variables import event_variables
from spanet_inference import ZHbbSpanetProcessor


class ZHbbSpanetTrainingProcessor(ZHbbSpanetProcessor):
    '''
    Same objects as the analysis (ZHbbSpanetProcessor: DeepJet AK4, ParticleNet XbbVsQCD and
    regressed mass, type-1 MET); after the preselection:
        ZH candidate variables (on all AK4 jets, as in the analysis selection), event-level SPANet
        inputs, optional AK4 cleaning (dR > 0.8 from the candidate) of the SPANet jets,
        process labels / gen flags (subsamples), prefilter to the training category,
        parton matching (MC, ../common/parton_matching.py)
    The event weight is the pocket-coffea weight (spanet_training.WEIGHTS); no object SFs.

    Workflow options (from spanet_training.py, on top of ZHbbSpanetProcessor's):
      clean_jets       - SPANet AK4 jets with dR > 0.8 from the candidate only
      matching         - {jet_dr, fatjet_dr, parton_copy}
      extra_b_samples  - samples whose extra hard-process b quarks get provenance 1 (4F ttbb)
      spanet_selection - spanet_training.SELECTION: the training category, applied as a prefilter
                         before the gen matching
      event_features   - {btag_wps, count_clean}: the event-level inputs (event_variables.py)
    '''

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.clean_jets = self.workflow_options.get("clean_jets", False)
        self.matching = self.workflow_options["matching"]
        self.extra_b_samples = self.workflow_options.get("extra_b_samples", [])
        self.event_feature_options = self.workflow_options["event_features"]
        selection = self.workflow_options.get("spanet_selection")
        self.prefilter_cut = get_spanet_selection(selection) if selection is not None else None

    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")

        # selection variables on all AK4 jets (as the analysis category)
        fj0 = self.events.FatJetSorted[:, 0]
        self.events['ZH_rap'] = vector.zip({"pt": fj0.pt, "eta": fj0.eta, "phi": fj0.phi, "mass": fj0.mass}).rapidity
        zh_helper(self.events)

        # event-level SPANet inputs, on all AK4 jets (the cleaning below only affects the jet inputs)
        opts = self.event_feature_options
        for name, values in event_variables(self.events.JetGood, ak.firsts(self.events.FatJetSorted),
                                            ak.firsts(self.events.LeptonGood), opts["btag_wps"],
                                            opts["count_clean"]).items():
            self.events[name] = values

        if self.clean_jets:
            dr = self.events.JetGood.delta_r(ak.firsts(self.events.FatJetSorted))
            self.events['JetGood'] = self.events.JetGood[(dr > 0.8) | ak.is_none(dr)]

        if self._isMC:
            match_gen_lep(self.events)
        self.label_process()          # process label + gen flags used by the subsamples

        if self.prefilter_cut is not None:
            self.apply_prefilter(self.prefilter_cut)

        if self._isMC:
            do_parton_matching(self.events, dict(self.matching,
                                                 extra_b_as_boson=self._sample in self.extra_b_samples))
