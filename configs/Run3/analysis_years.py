"""
Year-dependent configuration of the Run 3 analysis (2024 so far; add eras here for 2025).

  * eras  - data-taking periods as named in the datasets ('2024').
  * years - analysis years: one pickle directory and one datacard each (here identical to the eras).
"""
import re

# ---------------------------------------------------------------- eras
ERAS = ['2024']
ERA_ALIASES = {}
ERA_RE = '(?:' + '|'.join(sorted(ERAS + list(ERA_ALIASES), key=len, reverse=True)) + ')'


def find_era(s):
    """Era contained in a dataset name, e.g. 'TTtoLNu2Q_2024' -> '2024'."""
    m = re.search(ERA_RE, s)
    if m is None:
        raise ValueError(f"No era found in '{s}'")
    return ERA_ALIASES.get(m.group(), m.group())


# ---------------------------------------------------------------- analysis years
ERA_TO_YEAR = {e: e for e in ERAS}
ANALYSIS_YEARS = list(dict.fromkeys(ERA_TO_YEAR.values()))
YEARS = ['2024']           # years to process (one datacard per year)
YEAR_RE = '(?:' + '|'.join(sorted(ANALYSIS_YEARS, key=len, reverse=True)) + ')'


def find_year(s):
    m = re.search(YEAR_RE, s)
    if m is None:
        raise ValueError(f"No analysis year found in '{s}'")
    return m.group()


# ---------------------------------------------------------------- per-era settings
# Integrated luminosity in fb^-1 (golden JSON, runs C-I; same as pocket-coffea's lumi.yaml)
LUMI_FB = {'2024': 109.95}

# GloParT XbbVsQCD working points of the ZH candidate. Not defined centrally for Run 3 yet: to be
# derived from mistag rates (then e.g. {"2024": {"LP": ..., "MP": ..., "HP": ...}}).
XBB_WP_DICT = {}
XBB_WP = "LP"


def xbb_wp_cuts(era):
    """{'L', 'M', 'T'} GloParT XbbVsQCD thresholds of an era."""
    if era not in XBB_WP_DICT:
        raise KeyError(f"No GloParT Xbb working points defined for era '{era}' yet "
                       f"(XBB_WP_DICT in Run3/analysis_years.py)")
    wps = XBB_WP_DICT[era]
    return {"L": wps["LP"], "M": wps["MP"], "T": wps["HP"]}


def bbvl_wp(era):
    return xbb_wp_cuts(era)["L"]


# Pickle directory per analysis year
PICKLE_DIR = 'pickled/Nominal{year}'

# Luminosity uncertainties per analysis year: to be filled at the datacard stage (Run 3 LUM recommendation)
LUMI_UNCERTAINTIES = {}
