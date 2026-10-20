"""
pocket-coffea's JetsCalibrator with two kinds of variations that shift only some of the jets
(TOP recommendations for Run 2):

  AK4PFchs_JER_<region>Up/Down       JER varied only for jets in one (|eta|, pT) region of the
                                     JME Run 2 split (JER_REGIONS); replaces AK4PFchs_JERUp/Down
  AK4PFchs_JES_FlavorPure<F>Up/Down  pure-flavour JES source applied only to jets of that parton
                                     flavour (FLAVOUR_PARTONS); jets without partonFlavour (the
                                     low-pT CorrT1METJet of the type-1 MET) are not shifted

Every other jet keeps its nominal pT and mass. All other variations (JES sources, AK8 JER) are the
stock pocket-coffea ones. Regions and flavours use the nominally calibrated jet.
"""
import re

import awkward as ak
from pocket_coffea.lib.calibrators.common.common import JetsCalibrator

from analysis_years import JER_REGIONS

# pure-flavour JES source: |partonFlavour| values of the jets it applies to
FLAVOUR_PARTONS = {"Gluon": [21], "Quark": [1, 2, 3], "Charm": [4], "Bottom": [5]}

JER_SPLIT_RE = re.compile(r"^(AK4PFchs)_JER_(\w+?)(Up|Down)$")
FLAVOUR_RE = re.compile(r"^(AK4PFchs)_JES_FlavorPure(Gluon|Quark|Charm|Bottom)(Up|Down)$")


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
        variations = []
        for v in self._variations:
            m = re.match(r"^(AK4PFchs)_JER(Up|Down)$", v)
            if m:
                variations += [f"{m.group(1)}_JER_{region}{m.group(2)}" for region in JER_REGIONS]
            else:
                variations.append(v)
        self._variations = sorted(set(variations))

    def calibrate(self, events, orig_colls, variation, already_applied_calibrators=None):
        m_jer, m_flav = JER_SPLIT_RE.match(variation), FLAVOUR_RE.match(variation)
        if variation not in self._variations or not (m_jer or m_flav):
            return super().calibrate(events, orig_colls, variation, already_applied_calibrators)

        out = super().calibrate(events, orig_colls, "nominal", already_applied_calibrators)
        if m_jer:
            jet_type, region, direction = m_jer.groups()
            variation_type, select = "JER", (lambda jets: jer_region_mask(jets, region))
        else:
            jet_type, flavour, direction = m_flav.groups()
            variation_type, select = f"JES_FlavorPure{flavour}", (lambda jets: flavour_mask(jets, flavour))
        merged = self.jet_calib_param.get("merge_collections_for_variations", {}).get(self.year, {})
        for jt in merged.get(jet_type, [jet_type]):
            self.apply_masked_variation(out, jt, variation_type, direction.lower(), select)
        return out

    def apply_masked_variation(self, out, jet_type, variation_type, direction, select):
        """As JetsCalibrator.apply_variation, but only for the jets selected by `select`."""
        aliases = self.jet_calib_param.get("collection_name_alias", {}).get(self.year, {})
        jet_type_alias = aliases.get(jet_type, jet_type)
        jet_coll_name = self.jet_calib_param.collection[self.year][jet_type]
        if jet_coll_name is None:
            return
        jets = self.jets_calibrated[jet_coll_name]
        mask = select(jets)
        out[jet_coll_name]["pt"] = ak.where(mask, jets[f"pt_{variation_type}_{direction}"], jets.pt)
        out[jet_coll_name]["mass"] = ak.where(mask, jets[f"mass_{variation_type}_{direction}"], jets.mass)
        if self.jet_calib_param.sort_by_pt[self._year][jet_type_alias]:
            sorted_indices = ak.argsort(out[jet_coll_name]["pt"], axis=1, ascending=False)
            out[jet_coll_name] = out[jet_coll_name][sorted_indices]
            out[f"{jet_coll_name}_sortidx"] = sorted_indices
