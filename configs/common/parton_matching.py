"""
Gen matching for the SPANet training inputs (Run 2 and Run 3).

The partons are the hard-process quarks / gluons of GenPart and their provenance comes from the
ancestry (distinctParent), so it does not depend on the LHE record of the generator: the Run 2
functions (pocket_coffea.lib.parton_provenance) assume the Powheg ttH LHE ordering, which the
MadGraph TTZ-ZtoQQ-1Jets sample does not have. Provenance codes (as in Run 2):
  1 = H/Z -> qq (and, with extra_b_as_boson, extra hard-process b quarks: g -> bb in 4F ttbb)
  2 = b from the hadronic top    3 = b from the leptonic top
  4 = additional radiation       5 = quarks from the hadronic W      -1 = other
"""
import awkward as ak
import numpy as np

from pocket_coffea.lib.deltaR_matching import object_matching

BOSONS = (23, 25)


def assign_provenance(pdg, parent_pdg, grandparent_pdg, extra_b_as_boson=False):
    """Provenance of each parton from its pdgId and the pdgIds of its distinct parent and grandparent
    (jagged arrays, one list of partons per event)."""
    apdg, apar, agp = abs(pdg), abs(parent_pdg), abs(grandparent_pdg)
    is_quark = apdg <= 5
    from_boson = is_quark & ((apar == BOSONS[0]) | (apar == BOSONS[1]))
    from_w_top = is_quark & (apar == 24) & (agp == 6)
    b_from_top = (apdg == 5) & (apar == 6)
    # hadronic side = sign of the W decaying to quarks (W+ from t, W- from tbar); None if no hadronic W
    had_sign = ak.firsts(np.sign(parent_pdg[from_w_top]))
    had_sign = ak.fill_none(had_sign, 0)
    had_b = b_from_top & (np.sign(parent_pdg) == had_sign)
    extra = ((is_quark | (apdg == 21)) & ~from_boson & (apar != 6) & (apar != 24)
             & (apar != BOSONS[0]) & (apar != BOSONS[1]))

    prov = ak.values_astype(ak.full_like(pdg, -1), np.int64)
    prov = ak.where(extra, 4, prov)
    if extra_b_as_boson:
        prov = ak.where(extra & (apdg == 5), 1, prov)
    prov = ak.where(from_boson, 1, prov)
    prov = ak.where(had_b, 2, prov)
    prov = ak.where(b_from_top & ~had_b, 3, prov)
    prov = ak.where(from_w_top, 5, prov)
    return prov


def hard_partons(genpart, copy_flag="isFirstCopy", extra_b_as_boson=False):
    """Outgoing hard-process quarks and gluons with their provenance."""
    parent = genpart.distinctParent
    is_parton = (abs(genpart.pdgId) <= 5) | (genpart.pdgId == 21)
    # status 21 = incoming partons; an emission whose ancestors all have its pdgId (g <- g) has no
    # distinct parent and counts as additional radiation
    sel = genpart.hasFlags(["fromHardProcess", copy_flag]) & is_parton & (genpart.status != 21)
    partons = genpart[sel]
    parent = parent[sel]
    parent_pdg = ak.fill_none(parent.pdgId, 0)
    grandparent_pdg = ak.fill_none(parent.distinctParent.pdgId, 0)
    prov = assign_provenance(partons.pdgId, parent_pdg, grandparent_pdg, extra_b_as_boson)
    return ak.with_field(partons, prov, "provenance")


def hard_bosons(genpart):
    """Hard-process H / Z (last copy) decaying to quarks."""
    sel = ((abs(genpart.pdgId) == BOSONS[0]) | (abs(genpart.pdgId) == BOSONS[1])) \
        & genpart.hasFlags(["fromHardProcess", "isLastCopy"])
    bosons = genpart[sel]
    to_quarks = ak.any(abs(bosons.distinctChildren.pdgId) <= 5, axis=2)
    return bosons[to_quarks]


def match_fatjets(fatjets, bosons, dr_max):
    """matched / prov / dRMatchedJet of each fatjet w.r.t. the nearest H/Z (prov 1 if matched)."""
    dr = ak.min(fatjets.metric_table(bosons), axis=2)
    dr = ak.fill_none(dr, 999.)
    matched = dr < dr_max
    fatjets = ak.with_field(fatjets, matched, "matched")
    fatjets = ak.with_field(fatjets, ak.where(matched, 1, -1), "prov")
    return ak.with_field(fatjets, ak.where(matched, dr, -1.), "dRMatchedJet")


def do_parton_matching(events, options):
    """Sets Parton, PartonMatched, LeptonParton, BosonGen, nParton(Matched), nBosonPartonsInFatJet and
    the matched / prov / dRMatchedJet fields of JetGood and FatJetSorted.

    options: jet_dr, fatjet_dr, parton_copy, extra_b_as_boson
    """
    partons = hard_partons(events.GenPart, options["parton_copy"], options.get("extra_b_as_boson", False))
    matched_partons, matched_jets, dr = object_matching(partons, events.JetGood, dr_min=options["jet_dr"])

    jets = events.JetGood
    jets = ak.with_field(jets, ~ak.is_none(matched_partons, axis=1), "matched")
    jets = ak.with_field(jets, ak.fill_none(matched_partons.provenance, -1), "prov")
    jets = ak.with_field(jets, ak.fill_none(dr, -1.), "dRMatchedJet")
    events["JetGood"] = jets

    events["Parton"] = partons
    events["PartonMatched"] = matched_partons
    events["nParton"] = ak.num(partons, axis=1)
    events["nPartonMatched"] = ak.sum(jets.matched, axis=1)

    bosons = hard_bosons(events.GenPart)
    events["BosonGen"] = bosons
    events["FatJetSorted"] = match_fatjets(events.FatJetSorted, bosons, options["fatjet_dr"])

    # partons from the boson (prov 1) inside the candidate cone
    fj0 = ak.firsts(events.FatJetSorted)
    boson_partons = partons[partons.provenance == 1]
    events["nBosonPartonsInFatJet"] = ak.fill_none(ak.sum(boson_partons.delta_r(fj0) < 0.8, axis=1), 0)

    # outgoing LHE leptons / neutrinos
    lhe = events.LHEPart
    events["LeptonParton"] = lhe[(lhe.status == 1) & (abs(lhe.pdgId) > 10) & (abs(lhe.pdgId) < 17)]
