"""
AK8 jet mass scale and resolution (JMS, JMR) of the DNN-regressed fatjet mass (Run 2: ParticleNet,
Run 3: GloParT 2-prong), shared by the Run 2 and Run 3 analyses.

JMAR recommendation for DNN-regressed masses: no dedicated JMS or JMR correction (nominal SF = 1, the
regressed mass as is, in data and MC), and conservative uncertainties of 2% (JMS) and 3% (JMR). The
supporting study is Run 2; no Run 3 study exists, so the same values are used for Run 3.

  JMSUp / JMSDown   mass x 1.02 / 0.98
  JMRUp / JMRDown   resolution x 1.03 / 0.97 (apply_jms_jmr: hybrid smearing; the down variation
                    only narrows the gen-matched jets, unmatched jets cannot be unsmeared)

pocket-coffea runs one pass of the processor per shape variation that a calibrator declares. The
jet mass is not varied by a pocket-coffea calibrator: the processors apply JMS/JMR to the regressed
mass themselves (apply_object_preselection), reading the name of the current variation ("JMSUp", ...)
through jms_jmr_sfs. FatJetMassVariations only declares those four variations, so that pocket-coffea
runs them; it changes no collection. Select it with PC_JEC_VARIATION=JMSJMR.
"""
import zlib

import awkward as ak
import numpy as np
from pocket_coffea.lib.calibrators.calibrator import Calibrator
from pocket_coffea.lib.deltaR_matching import object_matching

FATJET_MASS_VARIATIONS = ["JMSUp", "JMSDown", "JMRUp", "JMRDown"]

# JMAR: 2% JMS, 3% JMR, no nominal correction
JMS_SF = {"JMSUp": 1.02, "JMSDown": 0.98}
JMR_SF = {"JMRUp": 1.03, "JMRDown": 0.97}

# Relative MC mass resolution of the stochastic smearing (unmatched jets), the baseline the JMAR
# team used in their derivations
SIGMA_MC = 0.10


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


def jms_jmr_sfs(variation):
    """(JMS, JMR) scale factors of a shape variation: 1 for the nominal and every other variation."""
    return JMS_SF.get(variation, 1.0), JMR_SF.get(variation, 1.0)


def chunk_rng(events):
    """Random generator seeded by the chunk (file + entry range): reproducible between runs
    and identical for every shape variation of the same chunk."""
    md = events.metadata
    key = f"{md.get('filename', '')}:{md.get('entrystart', 0)}:{md.get('entrystop', 0)}"
    return np.random.default_rng(zlib.crc32(key.encode()))


def apply_jms_jmr(events, fatjets, jms_sf, jmr_sf, sig_mc=SIGMA_MC, rng=None, mass=None):
    """
    JMS and JMR applied to the regressed mass `mass` of `fatjets` (default: fatjets.particleNet_mass).
    With jms_sf = jmr_sf = 1 the mass is returned unchanged.
    rng: numpy Generator for the stochastic smearing (see chunk_rng); global numpy RNG if None.
    """
    if mass is None:
        mass = fatjets.particleNet_mass
    if jms_sf == 1.0 and jmr_sf == 1.0:
        return mass

    # 1. Jet mass scale
    mass_jms = mass * jms_sf

    # 2. Match to GenJetAK8 (requirement: dR < R/2)
    matched_genjets, matched_fatjets, deltaR = object_matching(events.GenJetAK8, fatjets, dr_min=0.4)
    has_match = ~ak.is_none(matched_genjets, axis=1)
    gen_mass = ak.fill_none(matched_genjets.mass, 0.0)

    # 3. Deterministic smearing (matched)
    mass_deterministic = gen_mass + jmr_sf * (mass_jms - gen_mass)

    # 4. Stochastic smearing (unmatched), N(0, 1) per fatjet
    n_fatjets = int(ak.sum(ak.num(fatjets)))
    gauss = rng.standard_normal(n_fatjets) if rng is not None else np.random.normal(0, 1, n_fatjets)
    random_gauss = ak.unflatten(gauss, ak.num(fatjets))
    smear_factor = np.sqrt(np.maximum(jmr_sf ** 2 - 1.0, 0.0))
    mass_stochastic = mass_jms * (1.0 + random_gauss * sig_mc * smear_factor)

    # 5. Combine; no negative masses from extreme smearing
    mass_corr = ak.where(has_match, mass_deterministic, mass_stochastic)
    return ak.where(mass_corr < 0, 0.0, mass_corr)
