"""
Shape variations of the AK8 jet mass scale and resolution (JMS, JMR).

pocket-coffea runs one pass of the processor per shape variation that a calibrator declares.
The jet mass is not corrected by a pocket-coffea calibrator: the processor applies JMS/JMR to the
ParticleNet regressed mass itself (ZHbbSpanetProcessor.apply_object_preselection), and it reads the
name of the current variation ("JMSUp", ...) to pick the up/down correction.

This calibrator only declares those four variations, so that pocket-coffea runs them. It changes
no collection. Select it with PC_JEC_VARIATION=JMSJMR (inference_common.load_parameters).
"""
from pocket_coffea.lib.calibrators.calibrator import Calibrator

FATJET_MASS_VARIATIONS = ["JMSUp", "JMSDown", "JMRUp", "JMRDown"]


class FatJetMassVariations(Calibrator):
    name = "fatjet_mass"
    has_variations = True
    isMC_only = True
    calibrated_collections = []

    def __init__(self, params, metadata, do_variations=True, **kwargs):
        super().__init__(params, metadata, do_variations, **kwargs)
        self._variations = list(FATJET_MASS_VARIATIONS)

    def initialize(self, events):
        pass

    def calibrate(self, events, orig_colls, variation, already_applied_calibrators=None):
        return {}
