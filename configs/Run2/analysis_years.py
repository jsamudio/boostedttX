"""
Year-dependent configuration shared by coffea_to_pkl.py, make_datacard.py and datacard_shapes.py.

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


XBB_WP_OVERRIDE = 0.85      # e.g. 0.85 to test a looser cut (expected-sensitivity studies only)

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