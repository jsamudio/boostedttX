#                    #
##                  ##
######################
### Build datacard ###
######################
######################
#
# Usage:
#   python makeDatacard.py               # tt_B is a normal background
#   python makeDatacard.py --eft-ttbb    # tt_B declared as a signal so the EFT
#                                         # interferenceModel can scale it
# The --eft-ttbb cards get "eftTTB_" in their file names, so both versions can coexist.

import argparse
import numpy as np
import pandas as pd
import json
from pathos.multiprocessing import ProcessingPool as Pool
import re
import uproot
import os
import sys
from datacard_shapes import DataCardShapes, getZhbbWeight
from TH1 import export1d
from analysis_years import (YEARS, YEAR_RE, find_year, PICKLE_DIR, bbvl_wp,
                            LUMI_UNCERTAINTIES, jes_nuisances, JER_REGIONS)
from outvars import nominal_weight_factors, pdf_vars
from pdf_uncertainty import pdf_up_down, alphas_up_down, has_alphas

pt_bins = [0,200,300,450]
mass_bins = [50,80,105,145,200]
isblind = True

# =====================================================================
# Configuration
# =====================================================================
# Default for the tt_B switch (overridden by --eft-ttbb on the command line)
EFT_SCALE_TTBB = False

# Backgrounds in the card. Names must match the 'process' column / pickle file names.
# Add 'QCD' only if a QCD pickle is produced (the pickling script does not export it yet).
BACKGROUNDS = ['TTBar', 'tt_B', 'SingleTop', 'ttX', 'VJets']

# Normalization uncertainties of the minor backgrounds (lnN).
# PLACEHOLDERS -- set these to the values used in the previous analysis.
# (ttX and VJets have their own nuisances in add_Systematics: ttx_qsc, v_qsc, qqpdf)
BKG_XSEC_UNC = {'SingleTop': 1.10, 'QCD': 1.50}

# tt+2b (in tt_B) and tt+C (in TTBar) cross sections: events of the category scaled by these
# factors (up, down), as in the previous analysis. Categories: outvars.bkg_vars (tt_2b, tt_C).
TT_CATEGORY_XSEC_UNC = {'tt2bxsec': ('tt_2b', 'tt_B', 1.5, 0.5),
                        'ttCxsec':  ('tt_C', 'TTBar', 1.5, 0.5)}

# Theory event weights (scale, parton shower, PDF members) are clipped at these percentiles of
# each process's events in the pickle, against outlier weights (previous analysis)
THEORY_CLIP_PERCENTILES = (0.15, 99.85)
THEORY_CLIP_COLUMNS = [f'{s}_{ud}' for s in ['mu_r', 'mu_f', 'mu_rf', 'isr', 'fsr'] for ud in ['up', 'down']]

# Muon SFs split into a stat (per year) and a syst (correlated) nuisance each
MUON_SF_KEYS = ['id', 'iso', 'reco']
# Lepton SF extrapolation from the Z to the ttbar phase space (TOP): 0.5% upper limit, correlated
# across years, separate from the ID/iso SF nuisances
LEPTON_PHASESPACE_UNC = 1.005

# V+jets inclusive cross-section uncertainty: inclusive W -> l nu at NNLO (FEWZ 3.1),
# 20508.9 pb +165.7 -88.2 (scales) +-770.9 (PDF+alpha_s), CMS SM cross-section twiki (as before)
VJETS_SCALE_UNC = (1.008, 0.996)
VJETS_PDF_UNC   = 1.038

# Years, pickle paths, Xbb working points and luminosity uncertainties: see analysis_years.py

# Jet energy / mass shape uncertainties, from separately processed samples
# (run_pipeline.py --jec-variations JES_ALL JER JMSJMR -> pickles <process>_<tag>Up/Down.pkl).
USE_JET_VARIATIONS = False
# Processes the jet variations are applied to (as before); the others keep the nominal template
JET_SYST_PROCESSES = ['ttH', 'ttZ', 'TTBar', 'tt_B', 'SingleTop']


def jet_systematics(years):
    """(pickle tag, nuisance name) of every jet energy / mass nuisance. A name ending in the year
    is uncorrelated across years (written to that year's card only)."""
    out = jes_nuisances(years)                    # AK4 + AK8 shifted together (analysis_years.JES_SCHEME)
    for y in years:                               # resolution and mass: uncorrelated across years
        out += [(f'AK4PFchs_JER_{r}', f'jer_ak4_{r}_{y}') for r in JER_REGIONS]
        out += [('AK8PFPuppi_JER', f'jer_ak8_{y}'), ('JMS', f'jms_{y}'), ('JMR', f'jmr_{y}')]
    return out


# PDF / alpha_s shape uncertainties (requires the pdfweight_<i> columns in the pickles and the
# pdf/alphas factors in process_norms_<year>.json). The PDF set of a process follows from its number
# of members (pdf_uncertainty.py): 'pdf' for the 5FS Hessian set (signal, TTBar), 'pdf_ttbb' for the
# 4FS replica set (tt_B), uncorrelated as before; 'alphas' only for the Hessian set.
USE_PDF = True

# b-tag shape scale-factor sources. Columns: btag_sf<src>, btag_sf<src>_up, btag_sf<src>_down.
# Correlated sources get one nuisance in every year's card; the statistical ones one per year.
BTAG_SOURCES_CORRELATED   = ['hf', 'lf', 'cferr1', 'cferr2']
BTAG_SOURCES_UNCORRELATED = ['hfstats1', 'hfstats2', 'lfstats1', 'lfstats2']
# Denominator of the b-tag variations:
#   'source' -> btag_sf<src>_up / btag_sf<src>   (per-source central column)
#   'total'  -> btag_sf<src>_up / btag_sf        (if the per-source centrals are placeholders)
BTAG_DENOMINATOR = 'total'

def bbvl_threshold(df_, year):
    """Xbb cut per event from its era (a merged year like 2016 mixes eras with different
    working points); falls back to the year itself for pickles without an 'era' column."""
    if 'era' in df_:
        return df_['era'].map({e: bbvl_wp(e) for e in df_['era'].unique()})
    return bbvl_wp(year)

def cuts(df_, year):
    base = (
    (df_['n_ak4jets']   >= 5)       &
    (df_['n_b_outZH'] == 2) &
    (df_['ZH_bbvLscore'] >= bbvl_threshold(df_, year)) &
    (df_['MET_pt']      >= 20)            &
    (df_['ZH_M']        >= 50)            &
    (df_['ZH_M']        <= 200)
    )
    return base

class Systematic:
    '''
    Syntax: Systematic(Systematic Name,Systematic Type, Channel, Affected Processes, Value, Additional Information)
    '''
    dc_root_dir = ''
    datacard     = None
    allowed_processes = None
    all_years    = ['2018']

    def __init__(self, name, stype, process_ids, value, optvalue=None, info=None):
        self.name     = name
        self.stype    = stype
        self.ids      = process_ids
        # Year-specific nuisances carry the year in their name (written only to that year's card);
        # all others are written to every year's card (correlated across years).
        self.years     = re.findall(YEAR_RE, name) if name != 'jesHEMIssue' else ['2018']
        if len(self.years) == 0: self.years = list(self.all_years)

        if type(value) is list:
            self.value = {i:v for i,v in zip(self.ids,value)}
        elif type(value) is float or type(value) is int:
            self.value = value
        elif type(value) is dict:
            self.value = None
            self.years = list(value.keys())
        else:
            self.value = None

        self.optvalue = optvalue
        self.info     = '' if info is None else info

        if self.datacard is not None:
            for year in self.years:
                if type(value) is dict:
                    self.value = value[year]
                self.datacard[year].write(self.get_DC_line())

    @property
    def line(self):
        return '{0:14} {1:6}'.format(self.name,self.stype)

    @classmethod
    def set_dc_processes(cls,datacard, processes):
        cls.datacard = datacard
        cls.allowed_processes = processes
        cls.all_years = list(datacard.keys())
        # Process norms of the theory variations (written by make_process_norms.py)
        cls.p_norms = {}
        for y in cls.all_years:
            path = f'./process_norms/process_norms_{y}.json'
            if not os.path.exists(path):
                raise FileNotFoundError(f"{path} not found: run make_process_norms.py "
                                        f"(or run_pipeline.py norms) first")
            cls.p_norms.update(json.load(open(path, 'r')))

    def norm_factor(self, year, process, variation):
        """Factor of a theory variation from process_norms_<year>.json (tt_B is stored as 'ttbb')."""
        key = 'ttbb' if process == 'tt_B' else process
        try:
            return self.p_norms[year][key][variation]
        except KeyError:
            raise KeyError(f"No process norm for year={year}, process={key}, variation={variation} "
                           f"in process_norms_{year}.json: rerun make_process_norms.py") from None

    def get_DC_line(self):
        _line = self.line
        value = self.value
        optvalue = self.optvalue
        for p in self.allowed_processes:
            if p in self.ids:
                if type(self.value) is dict:
                    value = self.value[p]
                    if type(value) is list:
                        value, optvalue = value[0], value[1]
                    else:
                        value, optvalue = value, None
                if optvalue is None:
                    _line +='{0:12}'.format(str(value))
                else:
                    entry = '{1}/{0}'.format(str(value),str(optvalue))
                    _line += f'{entry:12}'
            else :
                _line +='{0:12}'.format('-')
        _line += '\t'+self.info+'\n'
        return _line

class ShapeSystematic(Systematic):
    bins_dict = None
    cut_op   = None
    df       = None

    def __init__(self, name, stype, subtype, process_ids, value, up=None, down=None, info=None, extraQC=False,
                 nominal=None, tag=None):
        super().__init__(name, stype, process_ids, value, info)
        self.up     = up
        self.down   = down
        self.tag    = tag or name       # 'jec': variation tag of the varied pickles
        self.subtype = subtype
        self.extraQC = extraQC
        self.nominal = nominal      # explicit nominal column; default: derived from the 'up' name

    def get_shape(self):
        fun_dict = {'mcstat' :self.makeMCStatHist,
                    'up/down'    :self.makeUpDownHist,
                    'normup/down':self.makeUpDownHist,
                    'ps'     :self.makeUpDownHist,
                    'pdfas'  :self.makeUpDownPdfasHist,
                    'erdOn'  :self.handleErdOn,
                    'jec'    :self.makeShapeVarHist}
        if self.subtype in fun_dict: fun_dict[self.subtype]()
        if self.extraQC and __name__ == '__main__': self.handle_extraQC()
        self.handleQC()
        return self.histos

    @classmethod
    def set_df_histos_histfuncs(cls, data_dict, hist_dict, sumw_sumw2):
        cls.data_dict = data_dict
        cls.histos = hist_dict
        cls.get_sumw_sumw2 = sumw_sumw2

    def makeShapeVarHist(self):
        """Up/down templates from separately processed (varied) samples, e.g. JES/JER/JMS/JMR."""
        for process in self.ids:
            for y in self.years:
                for ud in ['Up', 'Down']:
                    key = f'{process}_{y}_{self.tag}{ud}'
                    if key not in self.data_dict:
                        continue          # no events of this process: handleQC uses the nominal template
                    df = self.data_dict[key]
                    sumw, sumw2 = self.get_sumw_sumw2(df, getZhbbWeight(df, y), y)
                    self.histos[f'{process}_{y}_{self.name}{ud}'] = {'sumw': sumw, 'sumw2': sumw2}

    def handleErdOn(self):
        for process in self.ids:
            for y in self.years:
                df = self.data_dict[f'{process}_{y}_erdOn']
                weight =  getZhbbWeight(df,y)
                sumw, sumw2 = self.get_sumw_sumw2(df,weight,y)
                self.histos[f'{process}_{y}_{self.name}'] = {'sumw':sumw, 'sumw2': sumw2}

    @staticmethod
    def pdf_members(df):
        """Names of the pdfweight_<i> columns this process has (NaN columns = members it lacks)."""
        cols = [c for c in pdf_vars if c in df.columns and df[c].notna().any()]
        if any(df[c].isna().any() for c in cols):
            raise ValueError("Events of one process have different numbers of PDF members")
        return cols

    def makeUpDownPdfasHist(self):
        """PDF ('pdf_up'/'pdf_down') or alpha_s ('alphas_up'/'alphas_down') templates from the
        per-member templates (pdf_uncertainty.py), scaled by the process norms to remove the rate
        part (make_process_norms.py)."""
        kind = self.up.rsplit('_', 1)[0]                      # 'pdf' or 'alphas'
        for process in self.ids:
            for y in self.years:
                if f'{process}_{y}' not in self.data_dict: continue
                df = self.data_dict[f'{process}_{y}']
                nom_weight = getZhbbWeight(df,y)
                members = np.array([self.get_sumw_sumw2(df, nom_weight*df[c], y)[0] for c in self.pdf_members(df)])
                up, down = pdf_up_down(members) if kind == 'pdf' else alphas_up_down(members)
                up   = up   * self.norm_factor(y, process, self.up)
                down = down * self.norm_factor(y, process, self.down)
                self.histos[f'{process}_{y}_{self.name}Up']   = {'sumw': up,   'sumw2': up}
                self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw': down, 'sumw2': down}

    # Variations stored as standalone factors (no nominal column to divide out)
    STANDALONE_WEIGHTS = ('pdfWeight', 'mu_r', 'mu_f', 'mu_rf', 'isr', 'fsr') + tuple(TT_CATEGORY_XSEC_UNC)

    def makeUpDownHist(self):
        for process in self.ids:
            for y in self.years:
                if f'{process}_{y}' not in self.data_dict: continue
                df = self.data_dict[f'{process}_{y}']
                w_nom  = self.nominal or self.up.replace(re.search(r'(Up|up)\w*',self.up).group(),'').rstrip('_')
                # FIX: the old check ('mu_' in w_nom) also matched the muon scale factors
                # (mu_id_sf, mu_iso_sf, ...) and overwrote their nominal column with 1.0 in place,
                # which corrupted the nominal weight of every later systematic.
                standalone = (w_nom in self.STANDALONE_WEIGHTS) or (self.subtype == 'ps')
                weight =  getZhbbWeight(df,y)
                denom  = 1.0 if standalone else df[w_nom].replace(0, np.nan)
                for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
                    w_var = np.nan_to_num(weight * df[ud] / denom)
                    # 'normup/down' variations, and the PS variations of tt_B, are scaled by the
                    # process norm (make_process_norms.py). TTBar ISR/FSR ('ps') are not scaled.
                    if 'norm' in self.subtype or ('tt_' in process and self.subtype == 'ps'):
                        w_var = w_var * self.norm_factor(y, process, ud)
                    sumw, sumw2 = self.get_sumw_sumw2(df,w_var,y)
                    self.histos[f'{process}_{y}_{self.name}{ud_str}'] = {'sumw':sumw, 'sumw2':sumw2}

    def makeMCStatHist(self):
        for process in self.ids:
            for y in self.years:
                nom_hist = self.histos[f'{process}_{y}']
                stat_err = np.sqrt(nom_hist['sumw2'])
                mcstat_up, mcstat_down = nom_hist['sumw']+stat_err, nom_hist['sumw']-stat_err
                self.histos[f'{process}_{y}_{self.name}Up'] =   {'sumw' :mcstat_up, 'sumw2':mcstat_up}
                self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw' :mcstat_down, 'sumw2':mcstat_down}

    def handleQC(self):
        for process in self.ids:
            for y in self.years:
                nom_key = f'{process}_{y}'
                if nom_key not in self.histos:
                    continue

                nom     = np.float32(self.histos[nom_key]['sumw'])
                nom_err = np.sqrt(np.float32(self.histos[nom_key]['sumw2']))

                if re.search(r'(jms)|(jmr)', self.name) is not None and f'{nom_key}_{self.name}Up' not in self.histos:
                    tmp_name = re.search(r'(jms)|(jmr)', self.name).group()+f'_{y}'
                else:
                    tmp_name = self.name

                up_key = f'{nom_key}_{tmp_name}Up'
                down_key = f'{nom_key}_{tmp_name}Down'

                if up_key not in self.histos:
                    self.histos[up_key] = {'sumw': np.copy(nom), 'sumw2': np.copy(nom_err**2)}
                if down_key not in self.histos:
                    self.histos[down_key] = {'sumw': np.copy(nom), 'sumw2': np.copy(nom_err**2)}

                up   = np.float32(self.histos[up_key]['sumw'])
                down = np.float32(self.histos[down_key]['sumw'])

                up   = np.where(nom<nom_err, nom, up)
                down = np.where(nom<nom_err, nom, down)

                is_onesided = (( (up > nom) & (down > nom) ) | ( (up < nom) & (down < nom) ))
                geo_mean = np.sqrt( up * down )
                up , down = np.where(is_onesided, np.divide(up*nom, geo_mean) , up), np.where(is_onesided, np.divide(down*nom, geo_mean) , down)
                up, down = np.where((np.isinf(up)) | (np.isnan(up)), nom, up), np.where((np.isinf(down)) | (np.isnan(down)), nom, down)

                self.histos[up_key]['sumw']   = up
                self.histos[down_key]['sumw'] = down

    def handle_extraQC(self):
        for process in self.ids:
            for y in self.years:
                nom     = np.float32(self.histos[f'{process}_{y}']['sumw'])
                nom_err = np.sqrt(np.float32(self.histos[f'{process}_{y}']['sumw2']))
                if re.search(r'(jms)|(jmr)', self.name) is not None and f'{process}_{y}_{self.name}Up' not in self.histos:
                    tmp_name = re.search(r'(jms)|(jmr)', self.name).group()+f'_{y}'
                    self.histos[f'{process}_{y}_{self.name}Up'] = {'sumw':[],
                                                                   'sumw2':self.histos[f'{process}_{y}_{tmp_name}Up']['sumw2']}
                    self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw':[],
                                                                     'sumw2': self.histos[f'{process}_{y}_{tmp_name}Down']['sumw2']}
                else:
                    tmp_name = self.name
                up   = np.float32(self.histos[f'{process}_{y}_{tmp_name}Up']['sumw'])
                down = np.float32(self.histos[f'{process}_{y}_{tmp_name}Down']['sumw'])

                if 'hdamp' not in self.name and 'UE' not in self.name:
                    for i in range(nom.shape[-1]):
                        for j in range(nom.shape[0]):
                            if 'fsr' in self.name :
                                up[j,:,i]   = nom[j,:,i] * np.nansum(up[j,:,i])/np.nansum(nom[j,:,i])
                                down[j,:,i] = nom[j,:,i] * np.nansum(down[j,:,i])/np.nansum(nom[j,:,i])
                            else :
                                up[j,2:,i]   = nom[j,2:,i] * np.nansum(up[j,2:,i])/np.nansum(nom[j,2:,i])
                                down[j,2:,i] = nom[j,2:,i] * np.nansum(down[j,2:,i])/np.nansum(nom[j,2:,i])

                else:
                    for j in range(nom.shape[0]):
                        up[j,:,:]   = nom[j,:,:] * np.nansum(up[j,:,:])/np.nansum(nom[j,:,:])
                        down[j,:,:]   = nom[j,:,:] * np.nansum(down[j,:,:])/np.nansum(nom[j,:,:])

                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = up
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = down

class MakeDataCard:
    '''
    Handle creation and writing of datacard and shaped root file
    '''
    pt_bins = pt_bins
    mass_bins = mass_bins
    dc_dir = './datacards'

    def __init__(self, nn_name, hist2d, sig=['ttZ', 'ttH'], bkg=None, years=None, isblind=True,
                 eft_ttbb=EFT_SCALE_TTBB):
        years = list(YEARS) if years is None else list(years)
        bkg = list(BACKGROUNDS) if bkg is None else list(bkg)
        self.nn_name = nn_name
        self.hist2d = hist2d
        self.eft_ttbb = eft_ttbb
        # Tag includes the model name, and marks cards where tt_B is EFT-scaled
        self.tag = f'nominal_{self.nn_name}_' + ('eftTTB_' if eft_ttbb else '')
        self.sig_base = list(sig)
        self.bkg_base = list(bkg)
        self.sig = sig
        self.bkg = bkg
        self.data = ['data_obs']
        self.years = years
        self.isblind = isblind
        self.dc_bins = len(self.pt_bins[1:])
        self.cut = lambda x : x[self.nn_name] >= 0.0

        # Nominal weights and the up/down columns used for the weight systematics.
        # Missing columns are dropped silently by DataFrame.filter; systematics whose
        # columns are missing for a process are skipped for that process.
        weights = ['genWeight', 'norm_weight'] + nominal_weight_factors
        weight_sys = ['topptWeight_Up', 'topptWeight_Down',
                      'ele_id_sfup', 'ele_id_sfdown', 'ele_reco_sfup', 'ele_reco_sfdown',
                      'ele_trig_sfup', 'ele_trig_sfdown',
                      'mu_trig_sfup', 'mu_trig_sfdown',
                      'xbb_sfup', 'xbb_sfdown',
                      'puWeight_up', 'puWeight_down', 'l1prefire_sfup', 'l1prefire_sfdown',
                      'jet_puId_sfup', 'jet_puId_sfdown',
                      'isr_up', 'isr_down', 'fsr_up', 'fsr_down',
                      'mu_r_up', 'mu_r_down', 'mu_f_up', 'mu_f_down', 'mu_rf_up', 'mu_rf_down']
        weight_sys += [f'mu_{key}_sf{part}' for key in MUON_SF_KEYS
                       for part in ('statup', 'statdown', 'systup', 'systdown')]
        weight_sys += [f'btag_sf{src}{suffix}'
                       for src in BTAG_SOURCES_CORRELATED + BTAG_SOURCES_UNCORRELATED
                       for suffix in ('', '_up', '_down')]
        pdf_cols = list(pdf_vars) if USE_PDF else []

        self.weight_sys = weight_sys
        self.bkg_v = weights + weight_sys + pdf_cols + [self.nn_name, 'ZH_M', 'ZH_pt', 'process', 'sample',
                                                        'n_ak4jets', 'ZH_bbvLscore', 'n_b_outZH', 'MET_pt',
                                                        'era', 'year', 'tt_2b', 'tt_C']
        self.sig_v = self.bkg_v + ['genZHpt']

        self.accepted_sig  = [f'{s}{i}' for i in range(len(self.pt_bins)) for s in ['ttH','ttZ']]
        self.accepted_bkg  = list(bkg)
        self.accepted_data = ['data_obs']

        # Combine's interferenceModel only scales signal processes (index <= 0), so an
        # EFT-scaled background is declared as a signal in the card. This only changes
        # the process order/indices; templates, systematics and rateParams use names.
        self.eft_scaled_bkg = ['tt_B'] if eft_ttbb else []
        self.dc_sig   = self.accepted_sig + self.eft_scaled_bkg
        self.dc_bkg   = [b for b in self.accepted_bkg if b not in self.eft_scaled_bkg]
        self.dc_procs = self.dc_sig + self.dc_bkg

    # Moved inside the class to access self.nn_name and self.hist2d
    def get_sumw_sumw2(self, df, weights, year):
        ret_x_y = (lambda df: (df[self.nn_name].to_numpy(dtype='float'), df['ZH_M'].to_numpy(dtype='float')))
        sumw = np.array([
            self.hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]),
                                     weights=weights[df['pt_bin']==i_bin].astype('float'))[0]
            for i_bin in range(1,len(self.pt_bins))
        ], dtype = 'float')
        sumw2 = np.array([
            self.hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]),
                                     weights=np.power(weights[df['pt_bin']==i_bin],2).astype('float'))[0]
            for i_bin in range(1,len(self.pt_bins))
        ])
        return sumw, sumw2

    def makeDC(self):
        self.get_data()
        self.process_sig(self.data_dict)
        self.initialize_hists()
        self.initialize_roofile()
        self.initialize_datacard()
        self.setup_Systematics()
        self.add_Systematics()
        self.fill_roofile()
        self.close_roofile()
        self.close_dc()

    def get_data(self):
        self.data_dict = {}
        self.sig  = [f'{sig}_{y}'   for sig  in self.sig_base for y in self.years]
        self.bkg  = [f'{bkg}_{y}'   for bkg  in self.bkg_base for y in self.years]
        self.data = [f'{data}_{y}'  for data in ['data_obs'] for y in self.years]

        # Varied samples for the jet systematics (see jet_systematics)
        jet_tags = sorted({tag for tag, _ in jet_systematics(self.years)}) if USE_JET_VARIATIONS else []
        affected_mc = [p for p in self.sig_base + self.bkg_base if p in JET_SYST_PROCESSES]
        self.sys_samples = [f'{mc}_{y}_{tag}{ud}'
                            for mc in affected_mc
                            for y in self.years
                            for tag in jet_tags
                            for ud in ['Up', 'Down']]

        all_samples = self.sig + self.bkg + self.data + self.sys_samples

        pool = Pool(6)
        results = pool.map(self.worker, all_samples)
        pool.close()
        pool.join()
        pool.clear()

        for res in results:
            if res is None: continue
            for key, value in res.items():
                if key in self.data_dict:
                    self.data_dict[key] = pd.concat([self.data_dict[key],value], axis='rows', ignore_index=True)
                else:
                    self.data_dict[key] = value

        missing = [p for p in self.bkg_base if not any(k.startswith(f'{p}_') for k in self.data_dict)]
        if missing:
            print(f'[WARNING] no events found for backgrounds: {missing}')

    def worker(self, process):
        y = find_year(process)
        is_signal = any(sig_base in process for sig_base in ['ttH', 'ttZ'])
        p_vars = self.sig_v if is_signal else self.bkg_v
        if process in self.sys_samples:
            # varied samples only need the nominal weight (no up/down or PDF columns: saves memory)
            p_vars = [v for v in p_vars if v not in self.weight_sys and not v.startswith('pdfweight_')]
        if 'TTBar' in process : p_vars = p_vars + ['tt_type']
        if 'tt_B' in process : p_vars = p_vars + ['tt_type']
        if 'data'  in process : p_vars = p_vars + ['tt_type']

        return self.updatedict(process.replace(f'_{y}', '', 1), p_vars, y)   # tags may contain the year too

    def updatedict(self, p, v, y=''):
        filepath = f'{PICKLE_DIR.format(year=y)}/{p}.pkl'
        if not os.path.exists(filepath):
            return {}

        df = pd.read_pickle(filepath).filter(items=v)
        for col in [c for c in THEORY_CLIP_COLUMNS + list(pdf_vars) if c in df and df[c].notna().any()]:
            low, high = np.nanpercentile(df[col].to_numpy(dtype=float), THEORY_CLIP_PERCENTILES)
            df[col] = df[col].clip(low, high)
        # tt+2b / tt+C cross-section variations as event-weight factors (1 outside the category)
        for name, (flag, _, up, down) in TT_CATEGORY_XSEC_UNC.items():
            if flag in df:
                in_category = (df[flag] == True).to_numpy()
                df[f'{name}_up'] = np.where(in_category, up, 1.0)
                df[f'{name}_down'] = np.where(in_category, down, 1.0)
        df = df[cuts(df, y)]
        df = df.astype({k:'float32' for k in df if df[k].dtype == 'float64'})
        df.loc[:,'ZH_pt'] = df['ZH_pt'].clip(pt_bins[0]+1, pt_bins[-1]+1)
        df['pt_bin'] = pd.cut(df['ZH_pt'], bins=pt_bins+[np.inf], labels=[i_bin for i_bin in range(len(pt_bins))])

        df = df[self.cut(df)]
        group = df.groupby(by='process')

        sys_suffix = ""
        if "Up" in p or "Down" in p:
            base_proc = p.split('_')[0] if not 'tt_B' in p else 'tt_B'
            sys_suffix = "_" + p.replace(f"{base_proc}_", "")

        data_dict = {f"{n.replace('Data', 'data_obs')}_{y}{sys_suffix}": g for n,g in group}
        return data_dict

    def process_sig(self, data_dict_, y='2018'):
        sig_groups = ['ttZ', 'ttH']
        sig_names = []
        for g in sig_groups:
            sig_names = sig_names + re.findall(rf'{g}_{YEAR_RE}\w*', ' '.join(data_dict_.keys()))
        for sig_name in sig_names:
            s_pre = re.search(r'(ttH|ttZ)',sig_name).group()
            if sig_name not in data_dict_: continue
            df = data_dict_[sig_name]
            for i,_ in enumerate(self.pt_bins[:-1]):
                new_sig_name = f"{s_pre}{i}{sig_name.replace(s_pre, '')}"
                data_dict_[new_sig_name] = df[(df['genZHpt']>= self.pt_bins[i]) &
                                              (df['genZHpt'] < self.pt_bins[i+1])]
            data_dict_[f"{s_pre}{len(self.pt_bins) -1}{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            data_dict_.pop(sig_name)

    def initialize_hists(self):
        self.histos = {}
        for s,v in self.data_dict.items():
            print(s)
            y = find_year(s)
            w = getZhbbWeight(v,y) if 'data' not in s else np.ones_like((v[self.nn_name] if self.nn_name in v else v.iloc [:,0]).to_numpy())
            sumw, sumw2 = self.get_sumw_sumw2(v, w, y)
            self.histos[s] = {'sumw':sumw, 'sumw2':sumw2}

    def initialize_roofile(self):
        roo_dict = {}
        for y in self.years:
            roo_name = f'{self.dc_dir}/datacard_{self.tag}{y}.root'
            if os.path.exists(roo_name):
                os.system(f"rm {roo_name}")
            roo_dict[y] = uproot.create(roo_name)
        self.roo_dict = roo_dict

    def initialize_datacard(self):
        dc_dict = {y: open(f'{self.dc_dir}/datacard_{self.tag}{y}.txt', 'w') for y in self.years}
        for y,txt in dc_dict.items():
            txt.writelines([
                f'Datacard for {y}\n',
                'imax * number of bins\n',
                'jmax * number of processes minus 1\n',
                'kmax * number of nuisance paramerters\n',
                100*'-'+'\n',
                f'shapes data_obs * datacard_{self.tag}{y}.root $CHANNEL_data_obs\n',
                f'shapes * * datacard_{self.tag}{y}.root $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC\n',
                100*'-'+'\n',
                f"{'bin':20}{' '.join(['Zhpt'+str(i+1) for i in range(self.dc_bins)])}\n",
                f"{'observation':20}{self.dc_bins*'-1  '}\n",
                100*'-'+'\n',
                f"{'bin':20}{' '.join('Zhpt'+str(i+1) for i in range(self.dc_bins) for _ in self.dc_procs)}\n",
                f"{'process':20}{' '.join(s for _ in range(self.dc_bins) for s in self.dc_procs)}\n",
                f"{'process':20}{' '.join(str(i) for _ in range(self.dc_bins) for i in range(-len(self.dc_sig)+1, len(self.dc_bkg)+1))}\n",
                f"{'rate':20}{' '.join('-1' for _ in range(self.dc_bins) for _ in self.dc_procs)}\n",
                100*'-'+'\n'])
        self.dc_dict = dc_dict

    def write2dc(self, str2write):
        for y in self.years:
            self.dc_dict[y].write(str2write)

    def procs_with_columns(self, procs, cols, years=None):
        """Processes whose dataframes contain all `cols` in every requested year."""
        years = self.years if years is None else years
        ok = []
        for p in procs:
            dfs = [self.data_dict[f'{p}_{y}'] for y in years if f'{p}_{y}' in self.data_dict]
            if dfs and all(all(c in df.columns for c in cols) for df in dfs):
                ok.append(p)
        return ok

    def ttx_xsec_unc(self, year):
        try:
            unc = Systematic.p_norms[year]['ttX_xsec_unc']
        except KeyError:
            raise KeyError(f"No 'ttX_xsec_unc' in process_norms_{year}.json: rerun make_process_norms.py") from None
        if 'pdf_alphas' not in unc:
            raise KeyError(f"No ttX PDF uncertainty in process_norms_{year}.json: the outputs have no "
                           f"PDF weights (rerun the production, then make_process_norms.py)")
        return unc

    def pdf_member_counts(self, procs):
        """{process: number of PDF members} for the processes with PDF columns in every year."""
        counts = {}
        for p in procs:
            per_year = {len(ShapeSystematic.pdf_members(self.data_dict[f'{p}_{y}']))
                        for y in self.years if f'{p}_{y}' in self.data_dict}
            if per_year == {0} or not per_year:
                print(f'  [pdf] no PDF columns for {p} -> no PDF uncertainty')
                continue
            if len(per_year) > 1:
                raise ValueError(f'{p}: different numbers of PDF members across years: {per_year}')
            counts[p] = per_year.pop()
        return counts

    def add_Systematics(self):
        all_mc   = self.accepted_sig + self.accepted_bkg
        tth_sig  = [s for s in self.accepted_sig if 'ttH' in s]
        ttz_sig  = [s for s in self.accepted_sig if 'ttZ' in s]
        ttbar_mc = [p for p in ['TTBar', 'tt_B'] if p in self.accepted_bkg]
        minor_bkg = [p for p in self.accepted_bkg if p not in ttbar_mc]
        jec_mc   = [p for p in all_mc if re.sub(r'\d+$', '', p) in JET_SYST_PROCESSES]

        # Nuisance groups, written at the end (useful for stat/exp/theory breakdowns)
        groups = {'exp': [], 'theo': [], 'sig_rtheo': [], 'antisig_rtheo': [], 'ttbbnorm': []}

        def lnN(name, procs, value, optvalue=None, grps=()):
            procs = [p for p in procs if p in all_mc]
            if not procs:
                return
            Systematic(name, 'lnN', procs, value, optvalue)
            for g in grps:
                groups[g].append(name)

        def weight_shape(name, subtype, procs, up, down, grps, nominal=None, extraQC=False):
            cols = [up, down] + ([nominal] if nominal else [])
            procs_ok = self.procs_with_columns(procs, cols, re.findall(YEAR_RE, name) or None)
            skipped = [p for p in procs if p not in procs_ok]
            if skipped:
                print(f'  [{name}] missing {cols} for {skipped} -> not applied to them')
            if not procs_ok:
                return
            self.histos = ShapeSystematic(name, 'shape', subtype, procs_ok, 1, up, down, extraQC=extraQC,
                                          nominal=nominal).get_shape()
            for g in grps:
                groups[g].append(name)

        # ------------------------------------------------------------------
        self.write2dc('# Rate uncertainties\n')
        for lumi_name, per_year in LUMI_UNCERTAINTIES.items():
            per_year = {y: v for y, v in per_year.items() if y in self.years}
            if per_year:
                lnN(lumi_name, all_mc, per_year, grps=('exp',))
        lnN('lep_phasespace', all_mc, LEPTON_PHASESPACE_UNC, grps=('exp',))
        lnN('tth_ggpdf', tth_sig, 1.036,        grps=('sig_rtheo',))
        lnN('ttz_ggpdf', ttz_sig, 1.035,        grps=('sig_rtheo',))
        lnN('tth_qsc',   tth_sig, 1.058, 0.908, grps=('sig_rtheo',))
        lnN('ttz_qsc',   ttz_sig, 1.081, 0.907, grps=('sig_rtheo',))
        if 'recoeft' not in self.tag:
            lnN('tth_ggpdf0', ['ttH0'], 1.036,        grps=('antisig_rtheo',))
            lnN('ttz_ggpdf0', ['ttZ0'], 1.035,        grps=('antisig_rtheo',))
            lnN('tth_qsc0',   ['ttH0'], 1.058, 0.908, grps=('antisig_rtheo',))
            lnN('ttz_qsc0',   ['ttZ0'], 1.081, 0.907, grps=('antisig_rtheo',))
        lnN('ggpdf',  ttbar_mc, 1.042, grps=('theo',))
        lnN('tt_qsc', ttbar_mc, [[1.024, 0.965] for _ in ttbar_mc], grps=('theo',))
        lnN('qgpdf',  ['SingleTop'], 1.028, grps=('theo',))
        for p in minor_bkg:
            if p in BKG_XSEC_UNC:
                lnN(f'{p.lower()}_xsec', [p], BKG_XSEC_UNC[p], grps=('theo',))
        # ttX: from the samples' own scale / PDF weights (make_process_norms.py, 'ttX_xsec_unc');
        # VJets: inclusive W cross section. The PDF parts are qq-initiated, correlated (qqpdf).
        if 'ttX' in all_mc:
            ttx_unc = {y: self.ttx_xsec_unc(y) for y in self.years}
            lnN('ttx_qsc', ['ttX'], {y: {'ttX': [round(u['scale_up'], 3), round(u['scale_down'], 3)]}
                                     for y, u in ttx_unc.items()}, grps=('theo',))
        else:
            ttx_unc = {}
        lnN('v_qsc', ['VJets'], *VJETS_SCALE_UNC, grps=('theo',))
        qq_procs = [p for p in ['ttX', 'VJets'] if p in all_mc]
        lnN('qqpdf', qq_procs, {y: {**({'ttX': round(ttx_unc[y]['pdf_alphas'], 3)} if 'ttX' in qq_procs else {}),
                                    **({'VJets': VJETS_PDF_UNC} if 'VJets' in qq_procs else {})}
                                for y in self.years}, grps=('theo',))
        lnN('UE',         ['TTBar'], 1.01,  0.99,  grps=('theo',))
        lnN('hdamp',      ['TTBar'], 1.049, 0.954, grps=('theo',))
        lnN('hdamp_ttbb', ['tt_B'],  1.027, 0.974, grps=('theo',))

        self.write2dc('# Luminosity projection (freeze at fit time)\n')
        self.write2dc('lumiscale rateParam * * 1 [0,10]\n')

        # ------------------------------------------------------------------
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Shape uncertainties \n')

        # Jet energy / mass scale and resolution (separately processed samples)
        if USE_JET_VARIATIONS:
            for tag, name in jet_systematics(self.years):
                for y in re.findall(YEAR_RE, name) or self.years:
                    if not any(f'{p}_{y}_{tag}Up' in self.data_dict for p in jec_mc):
                        raise FileNotFoundError(f'No varied pickles for {tag} ({y}): run the {tag} variation '
                                                f'(run_pipeline.py --jec-variations) or set USE_JET_VARIATIONS = False')
                self.histos = ShapeSystematic(name, 'shape', 'jec', jec_mc, 1, tag=tag).get_shape()
                groups['exp'].append(name)

        # Experimental weight variations
        for y in self.years:
            for name, up, down, nominal in [
                (f'ele_id_sf_{y}',   'ele_id_sfup',   'ele_id_sfdown',   'ele_id_sf'),
                (f'ele_reco_sf_{y}', 'ele_reco_sfup', 'ele_reco_sfdown', 'ele_reco_sf'),
                (f'ele_trig_sf_{y}', 'ele_trig_sfup', 'ele_trig_sfdown', 'ele_trig_sf'),
                *[(f'mu_{key}_sf_stat_{y}', f'mu_{key}_sfstatup', f'mu_{key}_sfstatdown', f'mu_{key}_sf')
                  for key in MUON_SF_KEYS],
                (f'mu_trig_sf_{y}',  'mu_trig_sfup',  'mu_trig_sfdown',  'mu_trig_sf'),
                (f'xbb_sf_{y}',      'xbb_sfup',      'xbb_sfdown',      'xbb_sf'),
                (f'pu_{y}',          'puWeight_up',   'puWeight_down',   'puWeight'),
                (f'jet_puid_{y}',    'jet_puId_sfup', 'jet_puId_sfdown', 'jet_puId_sf'),
                (f'l1prefire_{y}',   'l1prefire_sfup', 'l1prefire_sfdown', 'l1prefire_sf'),
            ]:
                weight_shape(name, 'up/down', all_mc, up, down, ('exp',), nominal=nominal)

        # Muon reco / ID / iso SFs: systematic part correlated across years (MUO / TOP), the
        # statistical part per year (above)
        for key in MUON_SF_KEYS:
            weight_shape(f'mu_{key}_sf_syst', 'up/down', all_mc, f'mu_{key}_sfsystup', f'mu_{key}_sfsystdown',
                         ('exp',), nominal=f'mu_{key}_sf')

        # b-tag shape scale factors, per source
        btag_nom = (lambda src: 'btag_sf') if BTAG_DENOMINATOR == 'total' else (lambda src: f'btag_sf{src}')
        for src in BTAG_SOURCES_CORRELATED:
            weight_shape(f'btag_{src}', 'up/down', all_mc, f'btag_sf{src}_up', f'btag_sf{src}_down',
                         ('exp',), nominal=btag_nom(src))
        for y in self.years:
            for src in BTAG_SOURCES_UNCORRELATED:
                weight_shape(f'btag_{src}_{y}', 'up/down', all_mc, f'btag_sf{src}_up', f'btag_sf{src}_down',
                             ('exp',), nominal=btag_nom(src))

        # Modeling: top pT reweighting, parton shower, scales, PDF
        for name, (_, process, _, _) in TT_CATEGORY_XSEC_UNC.items():
            weight_shape(name, 'up/down', [process], f'{name}_up', f'{name}_down', ('theo',))
        weight_shape('toppt', 'up/down', ttbar_mc, 'topptWeight_Up', 'topptWeight_Down', ('theo',),
                     nominal='topptWeight')
        weight_shape('isr_tt',   'ps',          ['TTBar'], 'isr_up', 'isr_down', ('theo',), extraQC=True)
        weight_shape('fsr_tt',   'ps',          ['TTBar'], 'fsr_up', 'fsr_down', ('theo',), extraQC=True)
        weight_shape('isr_ttbb', 'ps',          ['tt_B'],  'isr_up', 'isr_down', ('theo',), extraQC=True)
        weight_shape('fsr_ttbb', 'ps',          ['tt_B'],  'fsr_up', 'fsr_down', ('theo',), extraQC=True)
        weight_shape('isr_tth',  'normup/down', tth_sig,   'isr_up', 'isr_down', ('theo',), extraQC=True)
        weight_shape('isr_ttz',  'normup/down', ttz_sig,   'isr_up', 'isr_down', ('theo',), extraQC=True)
        weight_shape('fsr_tth',  'normup/down', tth_sig,   'fsr_up', 'fsr_down', ('theo',), extraQC=True)
        weight_shape('fsr_ttz',  'normup/down', ttz_sig,   'fsr_up', 'fsr_down', ('theo',), extraQC=True)

        weight_shape('mu_r_tt',   'normup/down', ['TTBar'], 'mu_r_up', 'mu_r_down', ('theo',))
        weight_shape('mu_f_tt',   'normup/down', ['TTBar'], 'mu_f_up', 'mu_f_down', ('theo',))
        weight_shape('mu_r_tth',  'normup/down', tth_sig,   'mu_r_up', 'mu_r_down', ('theo',), extraQC=True)
        weight_shape('mu_f_tth',  'normup/down', tth_sig,   'mu_f_up', 'mu_f_down', ('theo',), extraQC=True)
        weight_shape('mu_r_ttz',  'normup/down', ttz_sig,   'mu_r_up', 'mu_r_down', ('theo',), extraQC=True)
        weight_shape('mu_f_ttz',  'normup/down', ttz_sig,   'mu_f_up', 'mu_f_down', ('theo',), extraQC=True)
        weight_shape('mu_r_ttbb', 'normup/down', ['tt_B'],  'mu_r_up', 'mu_r_down', ('theo',), extraQC=True)
        weight_shape('mu_f_ttbb', 'normup/down', ['tt_B'],  'mu_f_up', 'mu_f_down', ('theo',), extraQC=True)

        if USE_PDF:
            n_members = self.pdf_member_counts(self.accepted_sig + ttbar_mc)
            hessian  = [p for p, n in n_members.items() if has_alphas(n)]
            replicas = [p for p, n in n_members.items() if not has_alphas(n)]
            for name, procs, up, down in [('pdf',      hessian,  'pdf_up',    'pdf_down'),
                                          ('pdf_ttbb', replicas, 'pdf_up',    'pdf_down'),
                                          ('alphas',   hessian,  'alphas_up', 'alphas_down')]:
                if procs:
                    self.histos = ShapeSystematic(name, 'shape', 'pdfas', procs, 1, up, down).get_shape()
                    groups['theo'].append(name)

        # ------------------------------------------------------------------
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Float tt_B normalization\n')
        self.write2dc('CMS_ttbbnorm rateParam * tt_B 1 [0.0,5.0]\n')
        #groups['theo'].append('CMS_ttbbnorm')
        groups['ttbbnorm'].append('CMS_ttbbnorm')
        self.write2dc(100*'-'+'\n')
        self.write2dc('# MC Stats uncertainties\n')
        self.write2dc('* autoMCStats 10 0  1\n')
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Group definitions \n')
        for g, names in groups.items():
            if names:
                self.write2dc(f'{g} group = {" ".join(names)}\n')

    def setup_Systematics(self):
        # Same process order as the datacard header (dc_procs), so systematic columns line up
        process_line = np.array([self.dc_procs for _ in range(self.dc_bins)]).flatten()
        Systematic.set_dc_processes(self.dc_dict, process_line)
        ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos, sumw_sumw2=self.get_sumw_sumw2)

    @staticmethod
    def merge_last_mbin(pt_bin,a):
        a[pt_bin,:,-2] = a[pt_bin,:,-2] + a[pt_bin,:,-1]
        return a[pt_bin,:,:-1]

    @staticmethod
    def take_first_last_mbins(a):
        sb = np.stack([a[3:,0],a[3:,-1]], axis=-1)
        return np.append(a[0:3,:], sb).flatten()

    def fill_roofile(self):
        for p,v in self.histos.items():
            y = find_year(p)
            process = p.split(f'_{y}')[0]
            if process not in self.accepted_sig + self.accepted_bkg + self.accepted_data: continue
            if p.replace(f'_{y}', '') not in self.accepted_sig + self.accepted_bkg + self.accepted_data:
                sys     = p.replace(f'{process}_{y}', '')
            else:
                sys     = ''
            for pt_bin in range(v['sumw'].shape[0]):
                if pt_bin == 0:
                    if not self.isblind:
                        to_flat = (lambda a: self.take_first_last_mbins(self.merge_last_mbin(pt_bin,a)))
                    else:
                        to_flat = (lambda a: self.merge_last_mbin(pt_bin,a).flatten())
                else:
                    if not self.isblind:
                        to_flat = (lambda a : self.take_first_last_mbins(a[pt_bin,:,:]))
                    else:
                        to_flat = (lambda a : a[pt_bin,:,:].flatten())
                temp_dict = {'sumw' : to_flat(v['sumw'])}
                hist_name = f'Zhpt{pt_bin+1}_{process}{sys}'
                if 'sumw2' in v:
                    temp_dict['sumw2'] = to_flat(v['sumw2'])
                self.roo_dict[y][hist_name] = export1d(temp_dict, hist_name)

    def close_roofile(self):
        for roo in self.roo_dict:
            self.roo_dict[roo].close()

    def close_dc(self):
        for dc in self.dc_dict:
            self.dc_dict[dc].close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build the datacards.')
    parser.add_argument('--eft-ttbb', action='store_true', default=EFT_SCALE_TTBB,
                        help='Declare tt_B as a signal so the EFT interferenceModel can scale it '
                             '(cards get "eftTTB_" in their names)')
    args = parser.parse_args()

    from spanet_models import ACTIVE_MODELS     # one datacard per active SPANet model

    for model in ACTIVE_MODELS:
        nn_name = f"disc_{model}"
        print(f"\n{'='*50}")
        print(f" GENERATING DATACARDS FOR: {nn_name.upper()}"
              f"{'  (tt_B EFT-scaled, declared as signal)' if args.eft_ttbb else ''}")
        print(f"{'='*50}\n")

        hist2d = DataCardShapes(pt_bins, mass_bins, nn=nn_name, isblind=isblind, years=YEARS, common_binning=True) #False common binning for per-year
        hist2d.save_edges(f'nn_edges_{nn_name}.json')     # for NN_BINS in the EFT scaling
        from plot_nn_edges import plot_edges
        plot_edges(json.load(open(f'nn_edges_{nn_name}.json')), f'plots/nn_edges_{nn_name}.pdf')
        dc_maker = MakeDataCard(nn_name=nn_name, hist2d=hist2d, isblind=isblind, years=YEARS,
                                eft_ttbb=args.eft_ttbb)
        dc_maker.makeDC()
