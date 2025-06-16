'''
Contains commonly used functions for plotting
'''
import pandas as pd
import numpy as np
import os
import sys
import awkward as ak
from glob import glob

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter

def getZhbbBaseCuts(df_):
    base_cuts = (
        (df_['n_b_outZh']   >= 2)             &
        #(df_['n_lb_outZh']   >= 2)             &
        #(df_['n_b_outZh']   >= 1)             &
        (df_['n_ak4jets']   >= 5)             &
        #(df_['Zh_bbvLscore'] >= 0.8)          &
        (df_['Zh_bbvLscore'] >= 0.6)          &
        #Maybe unused( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        #Maybe unused(df_['passNotHadLep'] == 1) & # might add
        (df_['Zh_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['Zh_M']        >= 50)            &
        (df_['Zh_M']        <= 200))
    return base_cuts

def import_mpl_settings(i=1, width=1, length=1, disable_sansmath=False, no_figsize=False):
    import matplotlib.pyplot as plt
    plt.rc("font", size=10, family="sans-serif", **{"sans-serif" :
                                                    #[u'TeX Gyre Heros', u'Helvetica', u'Arial', u'Palatino']})
                                                    [u'TeX Gyre Heros', u'Helvetica', u'Arial']})
    plt.rc("xaxis", labellocation='right')
    plt.rc("yaxis", labellocation='top')
    plt.rc("legend", fontsize=10, scatterpoints=2, numpoints=1, borderpad=0.15, labelspacing=0.3,
           handlelength=0.7, handletextpad=0.25, handleheight=0.7, columnspacing=0.6,
           fancybox=False, edgecolor='none', borderaxespad=1)
    plt.rc("savefig", dpi=200)
    if not no_figsize :
        plt.rc("figure", figsize=(3.375*i*width, 3.375*(6./8.)*i*length), dpi=200)
    #plt.rc("mathtext", rm='sans')
    #if not disable_sansmath:
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble='\n'.join([
        r"\usepackage{amsmath}",
        #r"\usepackage{amssymb}",
        r"\usepackage{helvet}",
        #r"\usepackage{sansmath}",
        #r"\sansmath",
    ])+('\n'+r'\usepackage{sansmath}'+'\n'+r'\sansmath' if not disable_sansmath else ''))

class Plotter:
    '''
    master class containing main plotting functions
    grab files
    setup dependencies
    '''
    cut_func = staticmethod(getZhbbBaseCuts)
    fontsize = 10
    saveDir = 'pdf/'
    data_samples = None
    mc_samples = None
    year = None
    data_dict = None
    w_dict = None
    i_dict = None
    real_data = None
    HEM_opt = ''
    #

    def __init__(self, samples, kinem, bin_range,
                 xlabel=None, n_bins=20, bins=None,
                 doLog=True, doNorm=False, 
                 doShow = True, doSave = False,
                 doCuts=True,add_cuts=None,add_d_cuts=None,sepGenOpt=None,addData=False,addSoB=False,
                 alt_weight=None, show_int=False):
        
        import_mpl_settings(2, disable_sansmath=True)

        self.samples    = samples
        self.kinem      = kinem
        self.bin_range  = bin_range if bin_range else [bins[0],bins[-1]]
        self.xlabel     = kinem if xlabel is None else xlabel
        self.n_bins     = n_bins
        self.bins       = np.arange(bin_range[0],bin_range[-1]+((bin_range[-1]-bin_range[0])/n_bins) , (bin_range[-1]-bin_range[0])/n_bins) if bins is None else np.array(bins)
        self.bin_w      = np.linspace(bin_range[0],bin_range[-1],bins+1)[1:] - np.linspace(bin_range[0],bin_range[-1],bins+1)[:-1] if type(bins) is int  else self.bins[1:]-self.bins[:-1]
        self.doLog      = doLog
        self.doNorm     = doNorm
        self.doShow     = doShow
        self.doSave     = doSave
        self.addData    = addData
        self.addSoB     = addSoB
        self.doCuts     = doCuts
        self.add_cuts   = add_cuts if add_cuts is not None else ''
        self.add_d_cuts = add_d_cuts if add_d_cuts is not None else ''
        self.sepGenOpt  = sepGenOpt
        self.alt_weight = alt_weight
        self.show_int = show_int
        #
        self.prepData()

    def prepData(self):
            data = {sample: self.apply_cuts(self.data_dict[sample]) for sample in self.samples}
            if self.addData:
                data_dir = f'{self.pklDir}{self.year}/data_files/'
                self.real_data  = {'Data': np.clip(np.hstack(
                    [self.apply_data_cuts(sample,pd.read_pickle(f'{data_dir}{sample}_val.pkl'))[self.kinem].values 
                    for sample in self.data_samples]),self.bin_range[0],self.bin_range[-1])}
            #
            self.data = data
            
            self.sepData()
            #
            if self.alt_weight is None:
                default_weight = (lambda v_, obj: 
                                v_['weight']* np.sign(v_['genWeight'])
                                * v_['topptWeight']
                                * (v_['HEM_weight'] if obj.year+obj.HEM_opt == '2018' else 1.0)
                                * (v_['lep_trigeffsf'])
                                * v_['lep_sf']
                                #* v_['dak8md_bbvl_sf']
                                #* v_['BTagWeight'] 
                                #* v_['puWeight']  
                                #* (v_['PrefireWeight'] if obj.year != '2018' else 1.0))
                                )
            else:
                default_weight = self.alt_weight

            
            self.w_dict = {k: default_weight(v,self) for k,v in self.data.items()}


            self.i_dict     = {k: sum(v) for k,v in self.w_dict.items()}
            self.i_err_dict = {k: np.sqrt(sum(np.power(v,2))) for k,v in self.w_dict.items()}
            #
            self.data = {k: np.clip(v[self.kinem],self.bin_range[0],self.bin_range[-1]) for k,v in self.data.items()}

    @classmethod
    def load_data(cls,year='2017',HEM_opt='',samples=None,tag=None, addBSF=False, byprocess=False):
        cls.year = year
        cls.addBSF = addBSF
        cls.mc_samples = samples if samples is not None else cfg.All_MC
        cls.HEM_opt = HEM_opt
        cls.lumi = round(cfg.Lumi[year+HEM_opt],1)
        cls.data_dict = {}
        cls.tag = '' if tag is None else '_'+tag
        cls.byprocess = byprocess
        #
        pool = Pool()
        #
        if byprocess:
            pro_worker = partial(cls.pro_worker, pklDir=cls.pklDir, year=cls.year, tag=cls.tag)
            #results = pool.map(pro_worker, cls.mc_samples)
            results = map(pro_worker, cls.mc_samples)
        else:
            nom_worker = partial(cls.nom_worker, pklDir=cls.pklDir, year=cls.year, tag=cls.tag)
            results = pool.map(nom_worker, cls.mc_samples)
            #results = map(nom_worker, cls.mc_samples)
        for result in results:
            if result is None: continue
            for k,v in result.items():
                #print(k)
                if k in cls.data_dict:
                    cls.data_dict[k] = pd.concat([cls.data_dict[k],v], axis='rows', ignore_index=True)
                else:
                    cls.data_dict[k] = v
        #
        pool.close()
        
    @staticmethod
    def pro_worker(sample, pklDir=None, year=None, tag=None):
        try:
            df = pd.read_pickle(f'{pklDir}{year}/mc_files/{sample}{tag}_val.pkl')
        except:
            return None
        group = df.groupby(by='process')
        o_dict = {n:g for n,g in group}
        return o_dict
    @staticmethod
    def nom_worker(sample, pklDir=None, year=None, tag=None):
        try:
            df = pd.read_pickle(f'{pklDir}{year}/mc_files/{sample}{tag}_val.pkl')  
            return {sample:df}
        except:
            return None
    @classmethod
    def set_cut_func(cls, custom_cut_func):
        cls.cut_func = staticmethod(custom_cut_func)
    @classmethod
    def reset_cut_func(cls):
        cls.cut_func = staticmethod(getZhbbBaseCuts)


def get_sigbkg(filein, sig, bkg, vars_, genmatchreq, genmatch_sig=False):

    genweight_df = pd.DataFrame.from_dict(filein['sum_signOf_genweights'], orient='index')

    dfList = []

    #Make signal df

    for i in sig:
        tmp = []
        for j in filein['columns'][f'{i}'].keys():
            for var in vars_+[genmatchreq]:
                if (var == 'norm_weight'):
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                else:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                tmp.append(inner_list)
        tmp = np.transpose(np.asarray(tmp, dtype=object))
        tmpDF = pd.DataFrame(data=tmp, columns=vars_+[genmatchreq])
        dfList.append(tmpDF)
    s_df = pd.concat(dfList, ignore_index=True)

    dfList = []

    #Make bkg df

    for i in bkg:
        tmp = []
        for j in filein['columns'][f'{i}'].keys():
            for var in vars_+['tt_type']:
                if ((var == 'norm_weight') & ('TTbb' not in j)):
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                    inner_list = [i / genweight_df[0][f'{j}'] for i in inner_list]
                else:
                    inner_list = filein['columns'][f'{i}'][f'{j}']['btag_mask'][f'events_{var}'].value.tolist()
                tmp.append(inner_list)
        tmp = np.transpose(np.asarray(tmp, dtype=object))
        tmpDF = pd.DataFrame(data=tmp, columns=vars_+['tt_type'])
        dfList.append(tmpDF)
    b_df = pd.concat(dfList, ignore_index=True)

    b_df = b_df[(b_df['process'] == 'TTBar') | (b_df['process'] == 'tt_B')]

    '''
    Background is automatically trimmed to just TTBar and ttbb, and if genmatched
    signal is needed this can be done after this function.
    '''

    return s_df, b_df

def beginPlt(addSoB=False):
    if addSoB:
        fig, (ax, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})

    else:
        fig, ax = plt.subplots()
    fig.subplots_adjust(
        top=0.88,
        bottom=0.11,
        left=0.11,
        right=0.88,
        hspace=0.0 if addSoB else 0.2,
        wspace=0.2
    )

    if addSoB:
        return fig, ax, ax2

    else:
        return fig, ax

def endPlt(fig, ax, bin_range, doLog=False, doNorm=False):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    #CMSlabel(self.fig, self.ax, lumi=self.lumi)
    #plt.xlabel(self.xlabel+self.tag, fontsize = self.fontsize)
    #ax.set_ylabel(f"{'fraction of yield' if self.doNorm else 'Events'} / {(np.round(self.bin_w[0],2) if len(set(self.bin_w)) == 1 else 'bin')}")#fontsize = self.fontsize)

    ax.set_xlim(bin_range)
    if doLog: ax.set_yscale('log')
    handles, labels = ax.get_legend_handles_labels()
    #if dostat:
    #    hatch_patch = Patch(hatch=10*'X', label='Stat Unc.',  fc='w')
    #    handles = handles + [hatch_patch]
    #    labels  = labels + ['Stat Unc.']
    ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize=8)
    ax.set_ylim(ymin=(0 if not doLog else (.1 if not doNorm else .001)), # .1 -> .001
                     ymax=(ax.get_ylim()[1]*(30 if doLog else 1.50) if not doNorm else 1.05))
    #if self.doSave: plt.savefig(f'{self.saveDir}{self.xlabel}_.pdf', dpi = 300)
    #if self.doShow:
    #    plt.show()
    #    plt.close(self.fig)
    #plt.close('all') # for lpc

def formatSoBplot(ax2):
    ax2.errorbar(x=bin_c, y = sig_hist/np.sqrt(bkg_hist),
                      xerr=(edges[1:]-edges[:-1])/2,
                      fmt='.', color='k', label=r'S/$\sqrt{\mathrm{B}}$')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(which='both', direction='in', top=True, right=True)
    ax2.yaxis.set_label_coords(-0.07,0.35)
    ax2.set_ylabel(r'$\mathrm{S/}\sqrt{\mathrm{B}}$')
    ax2.grid(True)
