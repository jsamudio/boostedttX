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
# Integrated luminosity in fb^-1 (golden JSON, runs C-I); the weights take it from params/lumi.yaml
LUMI_FB = {'2024': 109.82}

# AK4 b-tag working points (UParT) whose information the analysis uses: the b-jet counts n_bjets_<WP>
# (SPANet) and the M requirement of the btag_mask category. The b-tag SFs (multi-WP method 1a,
# custom_scale_factors.sf_btag_multiwp) and the efficiency maps cover exactly these, loosest first.
BTAG_SF_WPS = ["L", "M", "T"]

# GloParT XbbVsQCD working points of the ZH candidate. Not defined centrally for Run 3. PROVISIONAL:
# LP from matching the Run 2 LP tt+LF / tt+C composition (scan_xbb_wp.py composition: 0.247 / 0.226),
# M / T as the xbb_wps of the SPANet model run3_refactor (spanet_models.py). The LP is the analysis
# cut in makeDatacard.cuts (bbvl_wp); the processor applies no Xbb cut in Run 3.
XBB_WP_DICT = {"2024": {"LP": 0.25, "MP": 0.877, "HP": 0.975}}
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

# Luminosity uncertainties per analysis year (lnN), {nuisance: {year: value}}, as in Run 2; a nuisance
# is correlated across the years it has an entry for. TO BE FILLED (LUM POG recommendation for 2024);
# makeDatacard.py warns while it is empty, e.g. {'lumi_13p6TeV_2024': {'2024': 1.0XX}}
LUMI_UNCERTAINTIES = {}
