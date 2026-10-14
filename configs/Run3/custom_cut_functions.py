"""Cuts of the Run 3 analysis (PUPPI MET, per-primary-dataset single-lepton triggers)."""

import awkward as ak
import numpy as np
from pocket_coffea.lib.cut_definition import Cut


def _no_none(mask, fill=False):
    return ak.where(ak.is_none(mask), fill, mask)


# ---------------------------------------------------------------- skim
def _precut(events, params, year, sample, **kwargs):
    """Loose skim on raw NanoAOD objects (before calibration and object selection)."""
    return _no_none((ak.num(events.Muon) + ak.num(events.Electron) >= 1)
                    & (ak.num(events.Jet) >= 3)
                    & (ak.num(events.FatJet) >= 1)
                    & (events.PuppiMET.pt >= 20))


precut = Cut(name="preSkimCut", params={}, function=_precut)


# ---------------------------------------------------------------- triggers
# Single-lepton triggers per era, split by the primary dataset (PD) that records them
HLT_PATHS = {
    "2024": {
        "egamma": ['HLT_Ele30_WPTight_Gsf', 'HLT_Ele115_CaloIdVT_GsfTrkIdT'],
        "muon":   ['HLT_IsoMu24', 'HLT_Mu50', 'HLT_HighPtTkMu100', 'HLT_CascadeMu100'],
    },
}
hlt_paths = {era: p["egamma"] + p["muon"] for era, p in HLT_PATHS.items()}

# Data samples -> trigger stream of their PD (EGamma0/1 -> "egamma", Muon0/1 -> "muon")
DATA_STREAMS = {"DATA_EGamma": "egamma", "DATA_Muon": "muon"}


def _hlt_or(events, paths):
    """OR of the given paths; paths missing in the NanoAOD count as False."""
    mask = np.zeros(len(events), dtype=bool)
    for path in paths:
        name = path.replace("HLT_", "")
        if name in events.HLT.fields:
            mask = mask | ak.to_numpy(events.HLT[name])
    return mask


def hlt_selection(events, params, year, sample=None, isMC=True, **kwargs):
    """Single-lepton trigger with per-PD logic, so no data event is counted twice:
      MC          : any e/gamma OR muon trigger
      DATA_Muon   : any muon trigger
      DATA_EGamma : any e/gamma trigger AND no muon trigger (those events are in the Muon PD)
    """
    if year not in HLT_PATHS:
        raise ValueError(f"No trigger paths defined for era '{year}'")
    egamma = _hlt_or(events, HLT_PATHS[year]["egamma"])
    muon = _hlt_or(events, HLT_PATHS[year]["muon"])
    if isMC:
        return egamma | muon
    stream = DATA_STREAMS.get(sample)
    if stream == "muon":
        return muon
    if stream == "egamma":
        return egamma & ~muon
    raise ValueError(f"Data sample '{sample}' has no trigger stream in DATA_STREAMS "
                     f"(known: {list(DATA_STREAMS)}); add it to avoid double counting")


hlt_cut = Cut(name="hlt_cut", params={}, function=hlt_selection)


# ---------------------------------------------------------------- event selections
def _event_selection(events, params, year, sample, **kwargs):
    """Analysis selection: 1 tight lepton, >= 5 AK4 jets, >= 1 AK8 jet, PUPPI MET > 20 GeV."""
    return _no_none((events.nJetGood >= 5)
                    & (events.nFatJetGood >= 1)
                    & (events.PuppiMET.pt >= 20)
                    & (events.nMuonGood + events.nElectronGood == 1))


event_selection = Cut(name="eventSelection", params={}, function=_event_selection)


def _trig_event_selection(events, params, year, sample, **kwargs):
    """Trigger-efficiency selection (orthogonal method): ttbar e-mu, exactly one tight electron
    and one tight muon, >= 4 AK4 jets with 1-2 b tags, >= 1 AK8 jet, PUPPI MET > 20 GeV."""
    return _no_none((events.nJetGood >= 4)
                    & (events.nFatJetGood >= 1)
                    & (events.PuppiMET.pt >= 20)
                    & (events.nMuonGood == 1)
                    & (events.nElectronGood == 1)
                    & (events.nbJetGood > 0)
                    & (events.nbJetGood <= 2))


trig_event_selection = Cut(name="trigEventSelection", params={}, function=_trig_event_selection)


def _btag_mask(events, params, year, sample, **kwargs):
    return _no_none(events.nbJetGood >= 2)


btag_mask = Cut(name="btagMask", params={}, function=_btag_mask)


def _spanet_selection(events, params, year, sample, **kwargs):
    """SPANet training category (spanet_training.SELECTION); None disables a requirement."""
    mask = np.ones(len(events), dtype=bool)
    if params.get("min_b_jets") is not None:
        mask = mask & (events.nbJetGood >= params["min_b_jets"])
    if params.get("min_ak4_jets") is not None:
        mask = mask & (events.n_ak4jets >= params["min_ak4_jets"])
    if params.get("n_b_outZH") is not None:
        mask = mask & (events.n_b_outZH == params["n_b_outZH"])
    if params.get("zh_mass") is not None:
        low, high = params["zh_mass"]
        mask = mask & (events.ZH_M >= low) & (events.ZH_M <= high)
    if params.get("min_xbb") is not None:
        mask = mask & (events.ZH_bbvLscore >= params["min_xbb"])
    return _no_none(mask)


def get_spanet_selection(selection):
    return Cut(name="spanetSelection", params=dict(selection), function=_spanet_selection)


def _veto_low_mass_pair(pair):
    """True for events WITHOUT an opposite-sign lepton pair of 0 < m < 12 GeV (J/psi, Upsilon)."""
    resonance = (pair.charge == 0) & (pair.mass < 12) & (pair.mass > 0)
    return ak.where(ak.is_none(resonance), False, ~resonance)


vetoMu = Cut(name="vetoMu", params={}, function=lambda events, **kw: _veto_low_mass_pair(events.mu_softmu))
vetoE = Cut(name="vetoE", params={}, function=lambda events, **kw: _veto_low_mass_pair(events.e_softe))


# ---------------------------------------------------------------- subsamples
def eq_genTtbarId_100(events, params, year, sample, **kwargs):
    """genTtbarId % 100 equal to params['genTtbarId'] (int) or any of them (iterable).
    0: tt+LF; 41-46: tt+c(c); 51-56: tt+b(b) (GenHFHadronMatcher)."""
    allowed_ids = [0, 41, 42, 43, 44, 45, 46, 51, 52, 53, 54, 55, 56]
    ids = params["genTtbarId"]
    ids = [ids] if isinstance(ids, int) else list(ids)
    if any(i not in allowed_ids for i in ids):
        raise ValueError(f"genTtbarId % 100 must be one of {allowed_ids}, got {ids}")
    mask = ak.zeros_like(events.event, dtype=bool)
    for i in ids:
        mask = mask | (events.genTtbarId % 100 == i)
    return mask


def get_genTtbarId_100_eq(genTtbarId, name=None):
    if name is None:
        ids = [genTtbarId] if isinstance(genTtbarId, int) else list(genTtbarId)
        name = "genTtbarId_100_eq_" + "_".join(str(i) for i in ids)
    return Cut(name=name, params={"genTtbarId": genTtbarId}, function=eq_genTtbarId_100)


def _event_flag(events, params, **kwargs):
    return _no_none(events[params["flag"]])


def get_flag_cut(flag):
    """Cut on a boolean event column set by the workflow (e.g. Zbb from genmatcher.match_gen_sig)."""
    return Cut(name=f"flag_{flag}", params={"flag": flag}, function=_event_flag)


genMatch = Cut(name="genMatch", params={},
               function=lambda events, **kw: _no_none(events["matchedGen_ZHbb_bb"]))
non_genMatch = Cut(name="non_genMatch", params={},
                   function=lambda events, **kw: _no_none(~events["matchedGen_ZHbb_bb"]))


# ---------------------------------------------------------------- helpers
def sortbyscore(coll, score):
    return coll[ak.argsort(coll[score], axis=1, ascending=False)]
