#                    #
##                  ##
######################
### Build datacard ###
######################
######################

import numpy as np
import pandas as pd
import json
from pathos.multiprocessing import ProcessingPool as Pool
import re
import uproot
import os
import sys
from datacard_shapes import DataCardShapes
from TH1 import export1d

pt_bins = [0,200,300,450]
#pt_bins = [200,300,450]
mass_bins = [50,80,105,145,200]
nn = 'nnscore'
isblind = True

def cuts(df_):
    base = (
    (df_['n_ak4jets']   >= 5)       &
    (df_['n_b_outZH'] == 2) &
    (df_['ZH_bbvLscore'] >= 0.9870) &
    (df_['ZH_pt']       >= 200)& # 200
    (df_['MET_pt']      >= 20)            &
    (df_['ZH_M']        >= 50)            &
    (df_['ZH_M']        <= 200)
    )
    return base

def get_sumw_sumw2(df, weights, year):
    ret_x_y = (lambda df: (df[nn].to_numpy(dtype='float'), df['ZH_M'].to_numpy(dtype='float')))
    sumw = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=weights[df['pt_bin']==i_bin].astype('float'))[0] 
        for i_bin in range(1,len(pt_bins))
    ], dtype = 'float')
    sumw2 = np.array([
        hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
                             weights=np.power(weights[df['pt_bin']==i_bin],2).astype('float'))[0]
        for i_bin in range(1,len(pt_bins))
    ])

    return sumw, sumw2


class MakeDataCard:
    '''
    Handle creation and writing of
    datacard, and the root file with
    shaped for template fit
    '''

    pt_bins = [0,200,300,450]
    #pt_bins = [200,300,450]
    mass_bins = [50,80,105,145,200]
    #mass_bins = [50,80,105,145]
    isblind = True
    weights = ['genWeight', 'norm_weight', 'topptWeight'] # add SFs here and then another list of the systematic weights
    #nn = ['ttzbb', 'tthbb', 'ttlf', 'ttbb', 'ttcc']
    nn = 'nnscore'
    bkg_v = weights + [nn,'ZH_M', 'ZH_pt', 'process', 'sample', 'n_ak4jets', 'ZH_bbvLscore', 'n_b_outZH', 'MET_pt']
    sig_v = bkg_v + ['genZHpt']
    accepted_sig  = [f'{s}{i}' for i in range(len(pt_bins)) for s in ['ttH','ttZ']]
    #accepted_sig = ['ttH', 'ttZ']
    accepted_bkg  = ['TTBar','tt_B'] # add others later
    accepted_data = ['data_obs']

    dc_dir = './datacards'
    tag = 'uscms_scaled3_'

    def __init__(self,
                 sig = ['ttZ', 'ttH'],
                 bkg = ['TTBar', 'tt_B'],
                 years = ['2017'],
                 cut = None,
                 isblind=True,
                 sumw_sumw2=get_sumw_sumw2):
        self.sig = sig
        self.bkg = bkg
        self.data = ['data_obs']
        self.years = years
        self.cut = (lambda x : x[self.nn] >= 0.0) if cut is None else cut
        self.isblind = isblind
        self.dc_bins = len(self.pt_bins[1:])
        self.get_sumw_sumw2 = sumw_sumw2

    def makeDC(self):
        self.get_data()
        self.process_sig(self.data_dict)
        self.initialize_hists()
        self.initialize_roofile()
        self.initialize_datacard()
        # add sys
        self.setup_Systematics()
        self.add_Systematics()
        self.fill_roofile()
        self.close_roofile()
        self.close_dc()
        

    def get_data(self):
        self.data_dict = {}
        self.sig  = [f'{sig}_{y}'   for sig  in self.sig  for y in self.years]
        self.bkg  = [f'{bkg}_{y}'   for bkg  in self.bkg  for y in self.years]
        self.data = [f'{data}_{y}'  for data in self.data for y in self.years]
        pool = Pool(6)
        all_samples = self.sig+self.bkg+self.data
        print(all_samples)
        results = pool.map(self.worker, all_samples)
        pool.close()
        for results in results:
            if results is None: continue
            for key, value in results.items(): # should be one
                print(key)
                if key in self.data_dict:
                    self.data_dict[key] = pd.concat([self.data_dict[key],value], axis='rows', ignore_index=True)
                else:
                    self.data_dict[key] = value
        del results, pool
        

    def worker(self, process):
        #y = re.search(r'201\d',process).group()
        y = '2017'
        p_vars = None
        p_vars = self.sig_v if process in self.sig else self.bkg_v
        p_vars = p_vars
        # add rate uncertanties to the processing
        if 'TTBar' in process : p_vars = p_vars + ['tt_type']
        if 'tt_B' in process : p_vars = p_vars + ['tt_type']
        if 'data'  in process : p_vars = p_vars + ['tt_type']
        #print(process.replace(f'_{y}', ''))
        return self.updatedict(process.replace(f'_{y}', ''), p_vars, y)

    def updatedict(self, p, v, y=''):
        if p in ['ttZ', 'ttH', 'tt_B', 'TTBar']:
            v += [''] # add in a bunch of the pdf weights
        #df = pd.read_pickle(f'pickled/SpanetInferenceTruncated_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/SpanetInferenceDoubleWithGenPt_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/SpanetInferenceAssignment_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/SpanetBalanced24SDM_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/SpanetNoTTCC_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/DNNInference_{p}.pkl').filter(items=v)
        #df = pd.read_pickle(f'pickled/XbbVsQCD_{p}.pkl').filter(items=v)
        df = pd.read_pickle(f'pickled/USCMSposter_inference_{p}.pkl').filter(items=v)
        print(p)
        df = df[cuts(df)]
        # Next add the mu_rf, isr/fsr, pdf 3 siga, etc. but we don't have it currently
        df = df.astype({k:'float32' for k in df if df[k].dtype == 'float64'}) # memory saving operation
        df.loc[:,'ZH_pt'] = df['ZH_pt'].clip(pt_bins[0]+1, pt_bins[-1]+1)
        df['pt_bin'] = pd.cut(df['ZH_pt'], bins=pt_bins+[np.inf], labels=[i_bin for i_bin in range(len(pt_bins))])
        #print(df)
        df = df[self.cut(df)]
        #df = df[(df[nn] > 0.0)])
        group = df.groupby(by='process')
        
        if 'data_obs' not in p: # Data -> data_obs
            for n,g in group:
                if n not in ['ttZ', 'ttH'] + self.accepted_bkg:
                    df.drop(g.index, inplace=True)
        group = df.groupby(by='process')
        del df
        # then extract systematics (for which there are none currently)
        sys = ''
        data_dict = {}
        # ... working on sys
        data_dict = {f"{n.replace('Data', 'data_obs')}_{y}{sys}": g for n,g in group}
        #data_dict = {f"{n}_{y}{sys}": g for n,g in group}
        print(data_dict)
        return data_dict

    def process_sig(self, data_dict_, y='2017'):
        sig_groups = ['ttZ', 'ttH']
        sig_names = []
        for g in sig_groups:
            sig_names = sig_names + re.findall(rf'{g}_201\d\w*', ' '.join(data_dict_.keys()))
        for sig_name in sig_names:
            s_pre = re.search(r'(ttH|ttZ)',sig_name).group()
            if sig_name not in data_dict_: continue
            print("made it")
            df = data_dict_[sig_name]
            for i,_ in enumerate(self.pt_bins[:-1]):
                new_sig_name = f"{s_pre}{i}{sig_name.replace(s_pre, '')}"
                data_dict_[new_sig_name] = df[(df['genZHpt']>= self.pt_bins[i]) &
                                              (df['genZHpt'] < self.pt_bins[i+1])]
            data_dict_[f"{s_pre}{len(self.pt_bins) -1}{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            #data_dict_[f"{s_pre}0{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            data_dict_.pop(sig_name)
        #print(data_dict_)
            
    def getZhbbWeight(self,df_, year):
        tot_weight = (((df_['norm_weight']/(41.529))*(137.596+62+120)) * np.sign(df_['genWeight']) * df_['topptWeight']) # and other weights
        return tot_weight

    def initialize_hists(self):
        self.histos = {}
        edges = []
        for s,v in self.data_dict.items():
            print(s)
            y = re.search(r'201\d',s).group()
            w = self.getZhbbWeight(v,y) if 'data' not in s else np.ones_like((v[nn] if nn in v else v.iloc [:,0]).to_numpy())
            sumw, sumw2 = self.get_sumw_sumw2(v, w, y)
            self.histos[s] = {'sumw':sumw,
                         'sumw2':sumw2}
        #print(self.histos)
    def initialize_roofile(self):
        roo_dict = {}
        for y in self.years:
            roo_name = f'{self.dc_dir}/datacard_{self.tag}{y}.root'
            if os.path.exists(roo_name):
                os.system(f"rm {roo_name}")
            roo_dict[y] = uproot.create(roo_name)
        self.roo_dict = roo_dict # self, that is

    def initialize_datacard(self):
        # 1 per year
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
                f"{'bin':20}{' '.join(['Zhpt'+str(i+1) for i in range(self.dc_bins) for _ in range(len(self.accepted_sig + self.accepted_bkg))])}\n",
                f"{'process':20}{' '.join(s for _ in range(self.dc_bins) for s in self.accepted_sig+self.accepted_bkg)}\n",
                f"{'process':20}{' '.join(str(i) for _ in range(self.dc_bins) for i in range(-len(self.accepted_sig)+1, len(self.accepted_bkg)+1))}\n",
                f"{'rate':20}{' '.join(str(-1) for _ in range(self.dc_bins) for s in self.accepted_sig + self.accepted_bkg)}\n",
                100*'-'+'\n'])
        #
        self.dc_dict = dc_dict

    def write2dc(self, str2write):
        for y in self.years:
            self.dc_dict[y].write(str2write)
    def add_Systematics(self):
        # SIG = ttZbb[0,1,2,3] ttHbb[0,1,2,3]
        # BKG = ttX, TTBar, old_tt_bb, VJets, other
        #process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        all_mc = self.accepted_sig + self.accepted_bkg
        all_but_ttbb = self.accepted_sig + ['TTBar'] #,'ttX','VJets','single_t']
        #tth_sig  = re.findall(r'ttH\d', ' '.join(self.accepted_sig))
        #ttz_sig  = re.findall(r'\w*ttZ\d', ' '.join(self.accepted_sig))
        tth_sig = [s for s in self.accepted_sig if 'ttH' in s]
        ttz_sig = [s for s in self.accepted_sig if 'ttZ' in s]
        ttbar_mc = ['TTBar','tt_B']
        #jec_mc   = ttbar_mc + ['single_t'] + tth_sig + ttz_sig
        jec_mc   = ttbar_mc + tth_sig + ttz_sig
        #
        #Systematic.set_dc_processes(self.dc_dict, process_line)
        self.write2dc(f'# Rate uncertainties\n')
        # new lumi
        #Systematic('lumi_13TeV_2016',       'lnN',  all_mc, 1.01)
        Systematic('lumi_13TeV_2017',       'lnN',  all_mc, 1.02)
        #Systematic('lumi_13TeV_2018',       'lnN',  all_mc, 1.015)
        #Systematic('lumi_13TeV_correlated', 'lnN',  all_mc, {'2016':1.006, '2017': 1.009, '2018': 1.02}) # 0.6, 0.9, 2.0
        #Systematic('lumi_13TeV_1718',       'lnN',  all_mc, {'2017':1.006, '2018': 1.002}) # 0.6, 0.2
        Systematic('lumi_13TeV_correlated', 'lnN',  all_mc, {'2017': 1.009}) # 0.6, 0.9, 2.0
        Systematic('lumi_13TeV_1718',       'lnN',  all_mc, {'2017':1.006}) # 0.6, 0.2
        # signal pdf qsc
        Systematic('tth_ggpdf', 'lnN', tth_sig, 1.036)       
        Systematic('ttz_ggpdf', 'lnN', ttz_sig, 1.035)       
        Systematic('tth_qsc' ,  'lnN', tth_sig, 1.058,0.908) 
        Systematic('ttz_qsc'  , 'lnN', ttz_sig, 1.081,0.907) 
        if 'recoeft' not in self.tag:
            Systematic('tth_ggpdf0', 'lnN', ['ttH0'], 1.036)       
            Systematic('ttz_ggpdf0', 'lnN', ['ttZ0'], 1.035)       
            Systematic('tth_qsc0' ,  'lnN', ['ttH0'], 1.058,0.908) 
            Systematic('ttz_qsc0'  , 'lnN', ['ttZ0'], 1.081,0.907) 
        # background pdf qcs
        Systematic('ggpdf', 'lnN', ttbar_mc, 1.042)          
        #Systematic('qqpdf', 'lnN', ['ttX','VJets'],
        #           [1.045, # ttX
        #            1.038, # VJets
        #        ])
        Systematic('qgpdf', 'lnN', ['single_t'], 1.028)
        # 
        Systematic('tt_qsc'   , 'lnN', ttbar_mc, [[1.024,0.965] for _ in ttbar_mc])
        #Systematic('ttx_qsc'  ,  'lnN', ['ttX'], 1.181, 0.875) 
        #Systematic('singlet_qsc'  ,  'lnN', ['single_t'], 1.031,0.979) 
        #Systematic('v_qsc'    , 'lnN', ['VJets'], 1.008, .996)#1.008, 0.996) # .821/1.24
        # Shape Systatics
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Shape uncertainties \n')
        #ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)#, self.hist3d, self.ptclip)
        # when defining shape, must incluse whether it is a mcsta, scale, or up/down syst
        Systematic('UE'         , 'lnN', ['TTBar'], 1.01, 0.99) # from .2 %
        Systematic('hdamp'      , 'lnN', ['TTBar'], 1.049, 0.954)
        Systematic('hdamp_ttbb' , 'lnN', ['tt_B'],  1.027, 0.974)
        #
        self.write2dc(100*'-'+'\n')
        #self.histos = ShapeSystematic('tt2bxsec', 'shape', 'up/down', ['tt_B'],  1, 'tt2bxsecWeight_Up', 'tt2bxsecWeight_Down').get_shape()
        #self.histos = ShapeSystematic('ttCxsec',  'shape', 'up/down', ['TTBar'], 1, 'ttCxsecWeight_Up',  'ttCxsecWeight_Down').get_shape()
        self.write2dc('# Float tt_B normalization\n') 
        self.write2dc('CMS_ttbbnorm rateParam * tt_B 1 [0.0,5.0]\n')
        self.write2dc(100*'-'+'\n')
        self.write2dc('# MC Stats uncertainties\n') 
        self.write2dc('* autoMCStats 10 0  1\n') 
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Group definitions \n') 
        self.write2dc('sig_rtheo group = tth_qsc ttz_qsc tth_ggpdf ttz_ggpdf\n')
        if 'recoeft' not in self.tag:
            self.write2dc('antisig_rtheo group = tth_qsc0 ttz_qsc0 tth_ggpdf0 ttz_ggpdf0\n')
        #self.write2dc('theo group = CMS_ttbbnorm tt2bxsec ttCxsec hdamp_ttbb hdamp UE toppt pdf_ttbb pdf alphas\n')
        self.write2dc('theo group = CMS_ttbbnorm hdamp_ttbb hdamp UE\n')
        #self.write2dc('theo group += mu_f_ttbb mu_r_ttbb mu_f_tt mu_r_tt mu_f_tth mu_r_tth mu_f_ttz mu_r_ttz\n')
        #self.write2dc('theo group += isr_ttbb fsr_ttbb isr_tt fsr_tt isr_tth fsr_tth isr_ttz fsr_ttz\n')
        #self.write2dc('theo group += tth_ggpdf ttz_ggpdf tth_qsc ttz_qsc ggpdf qqpdf qgpdf tt_qsc ttx_qsc singlet_qsc v_qsc\n')
        #
        
    def setup_Systematics(self):
        process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        Systematic.set_dc_processes(self.dc_dict, process_line)
        
    @staticmethod
    def merge_last_mbin(pt_bin,a):
        a[pt_bin,:,-2] = a[pt_bin,:,-2] + a[pt_bin,:,-1]
        return a[pt_bin,:,:-1]
        
    @staticmethod
    def take_first_last_mbins(a): # but not for the first 3 NN quantile
        sb = np.stack([a[3:,0],a[3:,-1]], axis=-1)
        return np.append(a[0:3,:], sb).flatten()
        
    def fill_roofile(self):
        for p,v in self.histos.items():
            print(p)
            y = re.findall(r'201\d',p)[0] # first instance of this should be the process year
            process = p.split(f'_{y}')[0]
            if process not in self.accepted_sig + self.accepted_bkg + self.accepted_data: continue
            if p.replace(f'_{y}', '') not in self.accepted_sig + self.accepted_bkg + self.accepted_data:
                # this should mean its a shape systematic or a process not in accepted
                sys     = p.replace(f'{process}_{y}', '') # should have format _sys
            else:
                sys     = ''
            #
            for pt_bin in range(v['sumw'].shape[0]):
                if pt_bin == 0: # 0,1,2 (-1)
                    if not self.isblind:
                        to_flat = (lambda a: self.take_first_last_mbins(self.merge_last_mbin(pt_bin,a)))
                    else:
                        to_flat = (lambda a: self.merge_last_mbin(pt_bin,a).flatten())
                else: # not the first pt bin
                    if not self.isblind:
                        to_flat = (lambda a : self.take_first_last_mbins(a[pt_bin,:,:]))
                    else:
                        to_flat = (lambda a : a[pt_bin,:,:].flatten())
                        #to_flat = (lambda a : self.merge_lastn_nnbin(pt_bin, a))
                temp_dict = {'sumw' : to_flat(v['sumw'])}#* (1 if y != '2017' else 3.3032)}#2.2967)} # to just scale to full run2
                #temp_dict = {'sumw' : to_flat(v['sumw'])* (1 if y != '2018' else cfg.Lumi['run2']/cfg.Lumi['2018'])} # to just scale to full run2
                hist_name = f'Zhpt{pt_bin+1}_{process}{sys}'
                if 'sumw2' in v:
                    temp_dict['sumw2'] = to_flat(v['sumw2']) 
                    #temp_dict['sumw2'] = to_flat(v['sumw2']) if 'tt_B' not in p else to_flat(v['sumw2'])/4 
                self.roo_dict[y][hist_name] = export1d(temp_dict, hist_name)
                
    def close_roofile(self):
        for roo in self.roo_dict:
            self.roo_dict[roo].close()
    
    def close_dc(self):
        for dc in self.dc_dict:
            self.dc_dict[dc].close()

class Systematic: # Class to handle Datacard systematics 
    ''' 
    Syntax: Systematic(Systematic Name,Systematic Type, Channel, Affected Processes, Value, Additional Information)
    '''
    #dc_root_dir = 'Higgs-Combine-Tool/'
    dc_root_dir = ''
    datacard     = None
    allowed_processes = None
    

    def __init__(self, name, stype, process_ids, value, optvalue=None, info=None):
        self.name     = name
        self.stype    = stype
        self.ids      = process_ids
        #self.channel  = channel
        self.years     = re.findall(r'201\d', name) if name != 'jesHEMIssue' else ['2018']
        if len(self.years) == 0: self.years = ['2017']
        # handle different value casts
        if type(value) is list:
            self.value = {i:v for i,v in zip(self.ids,value)}
        elif type(value) is float or type(value) is int:
            self.value = value
        elif type(value) is dict:
            self.value = None
            self.years = list(value.keys())
        else:
            self.value = None
        #self.value    = value if type(value) is not list else {i:v for i,v in zip(self.ids,value)}
        self.optvalue = optvalue
        self.info     = '' if info is None else info
        #
        if self.datacard is not None: 
            for year in self.years:
                if type(value) is dict:
                    self.value = value[year]
                self.datacard[year].write(self.get_DC_line()) # write to datacard file upon instance creation

    @property
    def line(self):
        return '{0:14} {1:6}'.format(self.name,self.stype)

    @classmethod
    def set_dc_processes(cls,datacard, processes):
        cls.datacard = datacard
        cls.allowed_processes = processes
        cls.p_norms = json.load(open(f'./process_norms/process_norms_ttbbw_run2.json','r'))

    def get_DC_line(self):
        _line = self.line
        value = self.value
        optvalue = self.optvalue
        for p in self.allowed_processes:
            #_process = p.replace('\t', '').replace(' ','' )# reformat process to exclude \t 
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
        
if __name__ == '__main__':
    hist2d = DataCardShapes(pt_bins,mass_bins,isblind=isblind)
    MakeDataCard(isblind=isblind).makeDC()


















    