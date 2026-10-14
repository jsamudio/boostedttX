"""
Type-1 PF MET for Run 2 UL, consistent between data and MC.

NanoAOD's MET is already type-1 corrected with the final UL JECs, in data and MC. pocket-coffea's
METCalibrator instead rebuilds MET from RawMET: in data (apply_data: false) it returns the RAW
MET, and in MC it sets the L1 (pileup offset) factor to 1, which is only right for PUPPI jets
(Run 3), not for the Run 2 CHS jets.

METType1PropagationCalibrator keeps NanoAOD's type-1 MET and only propagates the changes made by
the calibrators before it in the chain, i.e. the difference between the calibrated and the
original NanoAOD objects:
  * AK4 jets (JEC re-application, JER smearing, JES/JER variations), type-1 jet selection
    (corrected pT without muons > 15 GeV, EM fraction < 0.9); the L1 terms cancel in the difference
  * electron / muon scale corrections
In data the jets are not recalibrated (apply_jec_Data: false), so MET stays NanoAOD's type-1 MET
apart from the lepton scale corrections.

CALIBRATORS is the calibrator sequence to pass to the Configurator (default sequence with this
calibrator in place of METCalibrator).
"""
import awkward as ak
import numpy as np

from pocket_coffea.lib.calibrators.calibrator import Calibrator
from pocket_coffea.lib.calibrators.common.common import JetsCalibrator, ElectronsScaleCalibrator, MuonsCalibrator
from pocket_coffea.lib.jets import jets_in_original_order

T1_JET_PT_MIN = 15.0
T1_JET_EM_MAX = 0.9


def _sum_xy(pt, phi, mask):
    px = ak.sum(ak.where(mask, pt * np.cos(phi), 0.0), axis=1)
    py = ak.sum(ak.where(mask, pt * np.sin(phi), 0.0), axis=1)
    return px, py


def jet_t1_shift(jets_new, jets_old):
    """Change of the type-1 jet sum between two versions of the same AK4 jets (same order)."""
    nomu = 1.0 - jets_old.muonSubtrFactor
    em_ok = (jets_old.chEmEF + jets_old.neEmEF) < T1_JET_EM_MAX
    pt_new, pt_old = jets_new.pt * nomu, jets_old.pt * nomu
    nx, ny = _sum_xy(pt_new, jets_old.phi, em_ok & (pt_new > T1_JET_PT_MIN))
    ox, oy = _sum_xy(pt_old, jets_old.phi, em_ok & (pt_old > T1_JET_PT_MIN))
    return nx - ox, ny - oy


def lepton_shift(events, orig_colls, coll):
    key = f"{coll}.pt"
    if key not in orig_colls:
        return 0.0, 0.0
    lep = events[coll]
    dpt = lep.pt - orig_colls[key]
    return _sum_xy(dpt, lep.phi, ak.ones_like(dpt, dtype=bool))


class METType1PropagationCalibrator(Calibrator):
    name = "met_t1_propagation"
    has_variations = False
    isMC_only = False
    calibrated_collections = ["MET.pt", "MET.phi"]

    def initialize(self, events):
        pass

    def calibrate(self, events, orig_colls, variation, already_applied_calibrators=None):
        met = events["MET"]                      # NanoAOD type-1 MET (not touched before)
        mx = met.pt * np.cos(met.phi)
        my = met.pt * np.sin(met.phi)

        if "Jet" in orig_colls:                  # AK4 jets recalibrated / varied by JetsCalibrator
            dx, dy = jet_t1_shift(jets_in_original_order(events["Jet"]), orig_colls["Jet"])
            mx, my = mx - dx, my - dy
        for coll in ("Electron", "Muon"):        # lepton scale corrections
            dx, dy = lepton_shift(events, orig_colls, coll)
            mx, my = mx - dx, my - dy

        return {"MET.pt": np.hypot(mx, my), "MET.phi": np.arctan2(my, mx)}


# Default pocket-coffea sequence with the type-1 propagation in place of METCalibrator
CALIBRATORS = [JetsCalibrator, ElectronsScaleCalibrator, MuonsCalibrator, METType1PropagationCalibrator]
