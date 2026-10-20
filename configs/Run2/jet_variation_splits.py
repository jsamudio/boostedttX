"""
pocket-coffea's JetsCalibrator with two kinds of variations that shift only some of the jets
(TOP recommendations for Run 2):

  AK4PFchs_JER_<region>Up/Down       JER varied only for jets in one (|eta|, pT) region of the
                                     JME Run 2 split (JER_REGIONS); replaces AK4PFchs_JERUp/Down
  AK4PFchs_JES_FlavorPure<F>Up/Down  pure-flavour JES source applied only to jets of that parton
                                     flavour (FLAVOUR_PARTONS); jets without partonFlavour (the
                                     low-pT CorrT1METJet of the type-1 MET) are not shifted

Every other jet keeps its nominal pT and mass. The varied values come from pocket-coffea's own
apply_variation; all other variations (JES sources, AK8 JER) are the stock ones. Regions and
flavours use the nominally calibrated jet. Written for pocket-coffea's main branch (2026-10), where
the jets are re-sorted in JetsCalibrator.calibrate after apply_variation.
"""
import copy
import re

import awkward as ak
from pocket_coffea.lib.calibrators.common.common import JetsCalibrator

from analysis_years import JER_REGIONS

# pure-flavour JES source: |partonFlavour| values of the jets it applies to
FLAVOUR_PARTONS = {"Gluon": [21], "Quark": [1, 2, 3], "Charm": [4], "Bottom": [5]}

def split_selection(variation_type):
    """(stock variation, jet selection) for a split variation type, None for any other:
    'JER_<region>' -> ('JER', jets in the region); 'JES_FlavorPure<F>' -> (same, jets of flavour F)."""
    m = re.match(r"^JER_(\w+)$", variation_type)
    if m and m.group(1) in JER_REGIONS:
        return "JER", (lambda jets: jer_region_mask(jets, m.group(1)))
    m = re.match(r"^JES_FlavorPure(\w+)$", variation_type)
    if m and m.group(1) in FLAVOUR_PARTONS:
        return variation_type, (lambda jets: flavour_mask(jets, m.group(1)))
    return None


def jer_region_mask(jets, region):
    eta_min, eta_max, pt_min, pt_max = JER_REGIONS[region]
    abseta = abs(jets.eta)
    return (abseta >= eta_min) & (abseta < eta_max) & (jets.pt >= pt_min) & (jets.pt < pt_max)


def flavour_mask(jets, flavour):
    if "partonFlavour" not in jets.fields:
        return ak.zeros_like(jets.pt, dtype=bool)
    flav = abs(jets.partonFlavour)
    mask = ak.zeros_like(jets.pt, dtype=bool)
    for parton in FLAVOUR_PARTONS[flavour]:
        mask = mask | (flav == parton)
    return mask


class JetsCalibratorSplit(JetsCalibrator):
    # same name: it replaces the stock calibrator in the sequence (inference_common.CALIBRATORS)

    def initialize(self, events):
        super().initialize(events)
        # the stock AK4 JER variation is replaced by one variation per region
        variations = []
        for v in self._variations:
            m = re.match(r"^(AK4PFchs)_JER(Up|Down)$", v)
            if m:
                variations += [f"{m.group(1)}_JER_{region}{m.group(2)}" for region in JER_REGIONS]
            else:
                variations.append(v)
        self._variations = sorted(set(variations))

    def apply_variation(self, out, jet_type, variation_type, direction):
        """Stock apply_variation; for a split variation, the varied pt, mass and rawFactor are
        kept only for the selected jets (pocket-coffea re-sorts the jets afterwards, in calibrate)."""
        split = split_selection(variation_type)
        if split is None:
            return super().apply_variation(out, jet_type, variation_type, direction)
        stock_variation, select = split
        jet_coll_name = self.jet_calib_param.collection[self.year][jet_type]
        if jet_coll_name is None:
            return
        varied = {jet_coll_name: copy.copy(out[jet_coll_name])}
        super().apply_variation(varied, jet_type, stock_variation, direction)
        mask = select(self.jets_calibrated[jet_coll_name])
        nominal = out[jet_coll_name]
        for field in ("pt", "mass", "rawFactor"):
            out[jet_coll_name][field] = ak.where(mask, varied[jet_coll_name][field], nominal[field])
