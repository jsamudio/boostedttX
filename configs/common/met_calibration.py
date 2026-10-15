"""
Type-1 PF MET for Run 2 UL, recomputed from the recommended JECs (JME recipe).

JME: the type-1 MET stored in NanoAOD was computed with the JEC payload used at production and
must not be combined with AK4 jets corrected with a different payload. When the recommended JECs
are reapplied, the central type-1 MET and its JEC variations are recomputed from the same
corrected AK4 jets; only the AK4 jets entering the type-1 MET are propagated (not AK8).

METType1Calibrator therefore recomputes, for data and MC, nominal and every jet variation:

    MET = RawMET - sum_jets [ pT(L1L2L3[Res]) - pT(L1) ] * (1 - muonSubtrFactor)

over the AK4 type-1 jets, `Jet` and the low-pT `CorrT1METJet`, as corrected (and in MC smeared /
varied) by pocket-coffea's JetsCalibrator, with corrected pT (muon subtracted) > 15 GeV and, for
`Jet`, EM fraction < 0.9. The L1 factor comes from the same JEC payload (tag of the chunk: MC tag,
or the data tag of the run era). The lepton scale corrections of the calibrators before it are
propagated as well.

It replaces pocket-coffea's METCalibrator, which also starts from RawMET but uses an L1 factor of 1
(correct only for PUPPI jets) and returns the raw MET in data when apply_data is false.

Requires the recommended JECs to be reapplied to data as well (jets_calibration.apply_jec_Data),
so that `Jet` and `CorrT1METJet` carry the recommended payload in data and MC.

CALIBRATORS is the calibrator sequence to pass to the Configurator.
"""
import awkward as ak
import numpy as np

from pocket_coffea.lib.calibrators.calibrator import Calibrator
from pocket_coffea.lib.calibrators.common.common import JetsCalibrator, ElectronsScaleCalibrator, MuonsCalibrator
from pocket_coffea.lib.correction_cache import load_correction_set
from pocket_coffea.lib.jets import get_rho
from pocket_coffea.utils.utils import get_nano_version

T1_JET_PT_MIN = 15.0
T1_JET_EM_MAX = 0.9
T1_JET_TYPE = "AK4PFchs"                         # JEC payload of the type-1 jets (Run 2)
T1_COLLECTIONS = ("Jet", "CorrT1METJet")         # AK4 jets entering the type-1 MET


def _evaluate(corr, jets_flat, pt_raw_flat, rho_flat, run_flat):
    # "run": the Run 2 data JECs are one correction per year, binned in run number
    inputs = {"JetA": jets_flat.area, "JetEta": jets_flat.eta, "JetPt": pt_raw_flat, "Rho": rho_flat,
              "run": run_flat}
    return corr.evaluate(*[np.asarray(inputs[i.name], dtype=np.float64) for i in corr.inputs])


def _sum_xy(pt, phi, mask):
    px = ak.sum(ak.where(mask, pt * np.cos(phi), 0.0), axis=1)
    py = ak.sum(ak.where(mask, pt * np.sin(phi), 0.0), axis=1)
    return px, py


class METType1Calibrator(Calibrator):
    name = "met_type1_run2"
    has_variations = False
    isMC_only = False
    # Only existing fields: the CalibratorsManager cannot store originals of new fields
    calibrated_collections = ["MET.pt", "MET.phi"]

    def initialize(self, events):
        calib = self.params.jets_calibration.jet_types[T1_JET_TYPE][self.year]
        if self.isMC:
            tag = calib.jec_mc
        else:
            tag = calib.jec_data if isinstance(calib.jec_data, str) else calib.jec_data[self.metadata["era"]]
        cset = load_correction_set(calib.json_path)
        self.l1 = cset[f"{tag}_L1FastJet_{T1_JET_TYPE}"]
        full = f"{tag}_{calib.level}_{T1_JET_TYPE}"          # only for jets the JetsCalibrator did not touch
        self.l123 = cset.compound[full] if full in cset.compound.keys() else cset[full]
        self.rho = get_rho(events, get_nano_version(events, self.params, self.year))
        self.run = events.run

    def type1_sum(self, jets):
        """sum over type-1 jets of (pT_L1L2L3 - pT_L1) * (1 - muonSubtrFactor), as (px, py)."""
        counts = ak.num(jets)
        flat = ak.flatten(jets)
        rho = ak.flatten(ak.broadcast_arrays(self.rho, jets.eta)[0])
        run = ak.flatten(ak.broadcast_arrays(self.run, jets.eta)[0])
        if "pt" in jets.fields and "rawFactor" in jets.fields:   # corrected by the JetsCalibrator
            pt_corr = jets.pt
            pt_raw = jets.pt * (1.0 - jets.rawFactor)
        else:                                                     # e.g. CorrT1METJet left uncalibrated
            pt_raw = jets.rawPt
            pt_corr = pt_raw * ak.unflatten(_evaluate(self.l123, flat, ak.flatten(pt_raw), rho, run), counts)
        l1 = ak.unflatten(_evaluate(self.l1, flat, ak.flatten(pt_raw), rho, run), counts)
        nomu = 1.0 - jets.muonSubtrFactor
        pt_corr_nomu = pt_corr * nomu
        delta = pt_corr_nomu - pt_raw * l1 * nomu
        mask = pt_corr_nomu > T1_JET_PT_MIN
        if "chEmEF" in jets.fields and "neEmEF" in jets.fields:
            mask = mask & ((jets.chEmEF + jets.neEmEF) < T1_JET_EM_MAX)
        return _sum_xy(delta, jets.phi, mask)

    def calibrate(self, events, orig_colls, variation, already_applied_calibrators=None):
        raw = events["RawMET"]
        mx = raw.pt * np.cos(raw.phi)
        my = raw.pt * np.sin(raw.phi)
        for coll in T1_COLLECTIONS:
            if coll in events.fields:
                dx, dy = self.type1_sum(events[coll])
                mx, my = mx - dx, my - dy
        for coll in ("Electron", "Muon"):                        # lepton scale corrections
            key = f"{coll}.pt"
            if key in orig_colls:
                lep = events[coll]
                dx, dy = _sum_xy(lep.pt - orig_colls[key], lep.phi, ak.ones_like(lep.pt, dtype=bool))
                mx, my = mx - dx, my - dy
        return {"MET.pt": np.hypot(mx, my), "MET.phi": np.arctan2(my, mx)}


# Default pocket-coffea sequence with the Run 2 type-1 MET recomputation in place of METCalibrator
CALIBRATORS = [JetsCalibrator, ElectronsScaleCalibrator, MuonsCalibrator, METType1Calibrator]
