"""
Year-dependent configuration shared by the inference configs, outputToPkl.py, makeDatacard.py
and datacard_shapes.py.

Two levels:
  * eras  - data-taking periods as named in the datasets (2016_PreVFP, 2016_PostVFP, 2017, 2018).
            Luminosity and the Xbb working point are defined per era.
  * years - analysis years: one pickle directory and one datacard each.
            With COMBINE_2016 = True, both 2016 eras are merged into the analysis year '2016';
            every event keeps its era in the 'era' column, so era-specific cuts still apply.
"""
import re

# ---------------------------------------------------------------- eras
ERAS = ['2016_PreVFP', '2016_PostVFP', '2017', '2018']
ERA_ALIASES = {'2016APV': '2016_PreVFP', '2016': '2016_PostVFP'}   # dataset-name variants
ERA_RE = '(?:' + '|'.join(sorted(ERAS + list(ERA_ALIASES), key=len, reverse=True)) + ')'


def find_era(s):
    """Era contained in a dataset name, e.g. 'TTToHadronic_2016APV' -> '2016_PreVFP'."""
    m = re.search(ERA_RE, s)
    if m is None:
        raise ValueError(f"No era found in '{s}'")
    return ERA_ALIASES.get(m.group(), m.group())


# ---------------------------------------------------------------- analysis years
COMBINE_2016 = True
ERA_TO_YEAR = {e: ('2016' if (COMBINE_2016 and e.startswith('2016')) else e) for e in ERAS}
ANALYSIS_YEARS = list(dict.fromkeys(ERA_TO_YEAR.values()))          # e.g. ['2016', '2017', '2018']

# Years to process (one datacard per year)
YEARS = ['2018']           # e.g. ANALYSIS_YEARS for the full Run 2

YEAR_RE = '(?:' + '|'.join(sorted(ANALYSIS_YEARS, key=len, reverse=True)) + ')'


def find_year(s):
    """Analysis year contained in a process/histogram/nuisance name, e.g. 'ttZ0_2016_btag_sf_2016Up'."""
    m = re.search(YEAR_RE, s)
    if m is None:
        raise ValueError(f"No analysis year found in '{s}'")
    return m.group()


# ---------------------------------------------------------------- per-era settings
# Integrated luminosity per era in fb^-1 (Run 2 UL golden JSON)
LUMI_FB = {'2016_PreVFP': 19.52, '2016_PostVFP': 16.81, '2017': 41.48, '2018': 59.83}

# ZH candidate Xbb working points; the analysis uses the loose point (LP)
XBB_WP_DICT = {
    "2016_PreVFP":  {"LP": 0.9088, "MP": 0.9737, "HP": 0.9883},
    "2016_PostVFP": {"LP": 0.9137, "MP": 0.9735, "HP": 0.9883},
    "2017":         {"LP": 0.9105, "MP": 0.9714, "HP": 0.9870},
    "2018":         {"LP": 0.9172, "MP": 0.9734, "HP": 0.9880},
}
XBB_WP = "LP"


def xbb_wp_cuts(era):
    """{'L', 'M', 'T'} XbbVsQCD thresholds of an era: used for the xbb_wp_* SPANet inputs and the
    Xbb preselection in the processor (never affected by XBB_WP_OVERRIDE)."""
    wps = XBB_WP_DICT[era]
    return {"L": wps["LP"], "M": wps["MP"], "T": wps["HP"]}


XBB_WP_OVERRIDE = None      # e.g. 0.85 to test a looser cut (expected-sensitivity studies only)

def bbvl_wp(era):
    if XBB_WP_OVERRIDE is not None:
        return XBB_WP_OVERRIDE
    try:
        return XBB_WP_DICT[era][XBB_WP]
    except KeyError:
        raise KeyError(f"No Xbb {XBB_WP} working point for era '{era}'")


# Pickle directory per analysis year
PICKLE_DIR = 'pickled/Nominal{year}'

# ---------------------------------------------------------------- luminosity uncertainties
# Run 2 UL correlation scheme, per era. Verify against the current LUMI POG recommendation.
LUMI_UNCERTAINTIES_ERA = {
    'lumi_13TeV_2016':       {'2016_PreVFP': 1.010, '2016_PostVFP': 1.010},
    'lumi_13TeV_2017':       {'2017': 1.020},
    'lumi_13TeV_2018':       {'2018': 1.015},
    'lumi_13TeV_correlated': {'2016_PreVFP': 1.006, '2016_PostVFP': 1.006, '2017': 1.009, '2018': 1.020},
    'lumi_13TeV_1718':       {'2017': 1.006, '2018': 1.002},
}


def _lumi_per_year():
    out = {}
    for name, per_era in LUMI_UNCERTAINTIES_ERA.items():
        per_year = {}
        for era, val in per_era.items():
            y = ERA_TO_YEAR[era]
            if y in per_year and per_year[y] != val:
                raise ValueError(f"{name}: eras merged into {y} have different values")
            per_year[y] = val
        out[name] = per_year
    return out


LUMI_UNCERTAINTIES = _lumi_per_year()      # per analysis year, used by make_datacard.py


# ---------------------------------------------------------------- jet energy scale sources
# One shape variation per JES source (AK4 + AK8 jets shifted together), one production each
# (run_pipeline.py --jec-variations JES_ALL). The sources are checked against the JME files when
# the config is built (inference_common.py).
#   "full"     : the 27 individual sources (required by TOP)
#   "regrouped": the JME reduced set of 11 groups (6 correlated across years + 5 per year)
JES_SCHEME = "full"
# Flavour response (TOP default): FlavorQCD replaced by the pure-flavour sources, each applied to
# AK4 jets of that parton flavour only (jet_variation_splits.py). AK8 jets have no parton flavour
# in NanoAOD: they keep FlavorQCD as a separate nuisance (jes_FlavorQCD_AK8).
JES_FLAVOUR_SPLIT = True
JES_PURE_FLAVOURS = ["Gluon", "Quark", "Charm", "Bottom"]

JES_GROUPS = ["Absolute", "BBEC1", "EC2", "FlavorQCD", "HF", "RelativeBal"]       # correlated
JES_GROUPS_PER_YEAR = ["Absolute", "BBEC1", "EC2", "HF", "RelativeSample"]         # <group>_<year>

JES_SOURCES = [
    "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "FlavorQCD", "Fragmentation",
    "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpPtRef",
    "RelativeBal", "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
    "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeSample",
    "RelativeStatEC", "RelativeStatFSR", "RelativeStatHF", "SinglePionECAL", "SinglePionHCAL",
    "TimePtEta",
]
# Full sources uncorrelated across years (one nuisance per year); all others fully correlated.
# Same split as the regrouped scheme (the sources making up its <group>_<year> groups).
JES_SOURCES_UNCORRELATED = [
    "AbsoluteStat", "RelativeStatFSR", "TimePtEta",           # Absolute_<year>
    "RelativeJEREC1", "RelativePtEC1", "RelativeStatEC",       # BBEC1_<year>
    "RelativeJEREC2", "RelativePtEC2",                         # EC2_<year>
    "RelativeStatHF",                                          # HF_<year>
    "RelativeSample",                                          # RelativeSample_<year>
]


def jes_run_tags(scheme=JES_SCHEME):
    """Run tags (PC_JEC_VARIATION / --jec-variations) of every JES variation of a scheme.
    '_year' stands for the year of each era (e.g. JES_Regrouped_HF_year -> Regrouped_HF_2016)."""
    if scheme == "full":
        tags = [f"JES_{s}" for s in JES_SOURCES]
    else:
        tags = ([f"JES_Regrouped_{g}" for g in JES_GROUPS] +
                [f"JES_Regrouped_{g}_year" for g in JES_GROUPS_PER_YEAR])
    if JES_FLAVOUR_SPLIT:
        tags = [t for t in tags if not t.endswith("FlavorQCD")]
        tags += [f"JES_FlavorPure{f}" for f in JES_PURE_FLAVOURS] + ["JES_FlavorQCD_AK8"]
    return tags


def jes_source(tag, era):
    """JME source name of a run tag for one era, e.g. ('JES_Regrouped_HF_year', '2016_PreVFP') -> 'Regrouped_HF_2016'."""
    source = tag[len("JES_"):]
    if source == "FlavorQCD_AK8":
        return "FlavorQCD"
    return source[:-len("_year")] + "_" + era[:4] if source.endswith("_year") else source


def jes_jet_types(tag):
    """Jet types a JES run tag shifts: pure flavours only AK4 (+ type-1 MET jets), FlavorQCD_AK8
    only AK8, every other source AK4 and AK8 together (variation AK4AK8_JES_<source>)."""
    if tag.startswith("JES_FlavorPure"):
        return ["AK4PFchs"]
    if tag == "JES_FlavorQCD_AK8":
        return ["AK8PFPuppi"]
    return ["AK4PFchs", "AK8PFPuppi"]


def jes_nuisances(years, scheme=JES_SCHEME):
    """(variation tag in the pickles, nuisance name) of every JES nuisance. A name ending in the
    year is uncorrelated across years."""
    out = []
    if scheme == "full":
        out += [(f"AK4AK8_JES_{s}", f"jes_{s}_{y}") for s in JES_SOURCES_UNCORRELATED for y in years]
        out += [(f"AK4AK8_JES_{s}", f"jes_{s}") for s in JES_SOURCES if s not in JES_SOURCES_UNCORRELATED]
    else:
        out += [(f"AK4AK8_JES_Regrouped_{g}", f"jes_{g}") for g in JES_GROUPS]
        out += [(f"AK4AK8_JES_Regrouped_{g}_{y}", f"jes_{g}_{y}") for g in JES_GROUPS_PER_YEAR for y in years]
    if JES_FLAVOUR_SPLIT:
        out = [(t, n) for t, n in out if not t.endswith("FlavorQCD")]
        out += [(f"AK4PFchs_JES_FlavorPure{f}", f"jes_FlavorPure{f}") for f in JES_PURE_FLAVOURS]
        out += [("AK8PFPuppi_JES_FlavorQCD", "jes_FlavorQCD_AK8")]
    return out


# ---------------------------------------------------------------- jet energy resolution split
# AK4 JER: one variation per (|eta|, pT) region (JME Run 2 split; jet_variation_splits.py), one
# nuisance per region and year. name: (|eta| min, |eta| max, pT min, pT max)
JER_REGIONS = {
    "eta0to1p93":         (0.0, 1.93, 0.0, float("inf")),
    "eta1p93to2p5":       (1.93, 2.5, 0.0, float("inf")),
    "eta2p5to3pt0to50":   (2.5, 3.0, 0.0, 50.0),
    "eta2p5to3pt50toInf": (2.5, 3.0, 50.0, float("inf")),
    "eta3to5pt0to50":     (3.0, 5.0, 0.0, 50.0),
    "eta3to5pt50toInf":   (3.0, 5.0, 50.0, float("inf")),
}
