"""Workflow of the b-tag efficiency maps (Run 3): the AK4 jets of the analysis preselection, before the
b-tag requirement of the btag_mask category, in a histogram of hadron flavour x pT x |eta| x tag level
(BTAG_EFF_HISTOGRAMS in inference_common.py), for derive_btag_eff.py."""
import awkward as ak
import numpy as np

from pocket_coffea.utils.configurator import Configurator
from analysis_years import BTAG_SF_WPS
from spanet_inference import ZHbbRun3Processor


class ZHbbRun3BtagEffProcessor(ZHbbRun3Processor):
    '''
    Same objects and preselection as the analysis (no b-tag category). Adds to JetGood:
      abseta     - |eta|
      btag_level - number of the working points analysis_years.BTAG_SF_WPS the jet passes
                   (0 = untagged, len(BTAG_SF_WPS) = tagged at the tightest), from JetGood.btag_<WP>
    The ttbar subsamples (tt+LF / tt+C / tt+B) cut on genTtbarId directly; no other gen matching.
    '''

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)

    def process_extra_after_presel(self, variation):
        if not self._isMC:
            return
        jets = self.events.JetGood
        jets = ak.with_field(jets, np.abs(jets.eta), "abseta")
        level = sum(jets[f"btag_{wp}"] for wp in BTAG_SF_WPS)
        self.events['JetGood'] = ak.with_field(jets, level, "btag_level")
