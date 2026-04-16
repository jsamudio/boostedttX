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
from datacard_shapes import DataCardShapes, getZhbbWeight
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
    (df_['ZH_bbvLscore'] >= 0.9105) &
    #(df_['ZH_pt']       >= 200)& # 200
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
    weights = ['genWeight', 'norm_weight', 'topptWeight', 'ele_id_sf', 'ele_reco_sf', 'ele_id_sf', 'mu_id_sf', 'mu_iso_sf',
              'mu_reco_sf', 'bbtag_sf', 'btag_sf', 'puWeight'] # add SFs here and then another list of the systematic weights
    weight_sys = ['topptWeight_Up', 'topptWeight_Down', 'ele_id_sfup', 'ele_id_sfdown',
                #'ele_reco_sf', 
                'ele_reco_sfup', 'ele_reco_sfdown',
               #'ele_id_sf', 
                'ele_id_sfup', 'ele_id_sfdown',
               #'mu_id_sf', 
                  'mu_id_sfup', 'mu_id_sfdown',
               #'mu_iso_sf', 
                  'mu_iso_sfup', 'mu_iso_sfdown',
               #'mu_reco_sf', 
                  'mu_reco_sfup', 'mu_reco_sfdown',
               #'bbtag_sf', 
                  'bbtag_sfup', 'bbtag_sfdown',
               #'btag_sf',
               'btag_sfhf', 'btag_sfhf_up', 'btag_sfhf_down',
               'btag_sflf', 'btag_sflf_up', 'btag_sflf_down',
               'btag_sfhfstats1', 'btag_sfhfstats1_up', 'btag_sfhfstats1_down',
               'btag_sfhfstats2', 'btag_sfhfstats2_up', 'btag_sfhfstats2_down',
               'btag_sflfstats1', 'btag_sflfstats1_up', 'btag_sflfstats1_down',
               'btag_sflfstats2', 'btag_sflfstats2_up', 'btag_sflfstats2_down',
               'btag_sfcferr1', 'btag_sfcferr1_up', 'btag_sfcferr1_down',
               'btag_sfcferr2', 'btag_sfcferr2_up', 'btag_sfcferr2_down',
               #'puWeight', 
                'puWeight_up', 'puWeight_down',
               'isr_up', 'isr_down', 'fsr_up', 'fsr_down',
               'mu_r_up', 'mu_r_down',
               'mu_f_up', 'mu_f_down',
               'mu_rf_up', 'mu_rf_down',]
    #nn = ['ttzbb', 'tthbb', 'ttlf', 'ttbb', 'ttcc']
    nn = 'nnscore'
    bkg_v = weights + weight_sys + [nn,'ZH_M', 'ZH_pt', 'process', 'sample', 'n_ak4jets', 'ZH_bbvLscore', 'n_b_outZH', 'MET_pt']
    sig_v = bkg_v + ['genZHpt']
    accepted_sig  = [f'{s}{i}' for i in range(len(pt_bins)) for s in ['ttH','ttZ']]
    #accepted_sig = ['ttH', 'ttZ']
    #accepted_bkg  = ['TTBar','tt_B', 'VJets'] # add others later
    accepted_bkg  = ['TTBar','tt_B']#, 'VJets'] # add others later
    accepted_data = ['data_obs']
    all_sys_samples = ['ttZ_jesAbsoluteStatUp', 'ttZ_jesAbsoluteStatDown']

    dc_dir = './datacards'
    tag = 'run3test_'

    def __init__(self,
                 sig = ['ttZ', 'ttH'], #+['ttZ_jesAbsoluteStatUp', 'ttZ_jesAbsoluteStatDown'],
                 bkg = ['TTBar', 'tt_B'], #, 'VJets'],
                 years = ['2024'],
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
        #print(all_samples)
        results = pool.map(self.worker, all_samples)
        pool.close()
        for results in results:
            if results is None: continue
            for key, value in results.items(): # should be one
                #print(key)
                if key in self.data_dict:
                    self.data_dict[key] = pd.concat([self.data_dict[key],value], axis='rows', ignore_index=True)
                else:
                    self.data_dict[key] = value
        del results, pool
        

    def worker(self, process):
        #y = re.search(r'201\d',process).group()
        y = '2024'
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
        if p in self.all_sys_samples:
            v = [var for var in v if var not in self.weight_sys] # save memory
            df = pd.read_pickle(f'pickled/topitx/{p}_val.pkl').filter(items=v)
            print(df.keys())
        else:
            if p in ['ttZ', 'ttH', 'tt_B', 'TTBar']: #, 'VJets']:
                v += [''] # add in a bunch of the pdf weights
            #df = pd.read_pickle(f'pickled/SpanetInferenceTruncated_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/SpanetInferenceDoubleWithGenPt_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/SpanetInferenceAssignment_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/SpanetBalanced24SDM_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/SpanetNoTTCC_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/DNNInference_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/XbbVsQCD_{p}.pkl').filter(items=v)
            #df = pd.read_pickle(f'pickled/DataCardSys_{p}.pkl').filter(items=v)
            df = pd.read_pickle(f'pickled/Inference_{p}.pkl').filter(items=v)
            #print(p)
            #print(df)
            df = df[cuts(df)]
            #print(df)
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
       # del df #FIXME
        # then extract systematics (for which there are none currently)
        sys = ''
        data_dict = {}
        if p in self.all_sys_samples: 
            #print(p)
            if __name__ != '__main__': return data_dict
            sys =  '_'+p.split('_')[-1] # format [process_name]_[systype]
            #print(sys)
            if 'sys' in p: # hdamp and UE here
                data_dict.update({f"{n}_{y}_hdamp{'_ttbb' if 'ttbb' in p else ''}Up": g[g['sample'].str.contains('hdampUp')] for n,g in group})
                data_dict.update({f"{n}_{y}_hdamp{'_ttbb' if 'ttbb' in p else ''}Down": g[g['sample'].str.contains('hdampDown')] for n,g in group})
                if 'ttbb' not in p:
                    data_dict.update({f'{n}_{y}_UEUp': g[g['sample'].str.contains('UEUp')] for n,g in group})
                    data_dict.update({f'{n}_{y}_UEDown': g[g['sample'].str.contains('UEDown')] for n,g in group})
                return data_dict
        # ... working on sys
        data_dict = {f"{n.replace('Data', 'data_obs')}_{y}{sys}": g for n,g in group}
        #data_dict = {f"{n}_{y}{sys}": g for n,g in group}
        #print(data_dict)

        if not data_dict:
            process_name = p.replace('Data', 'data_obs')
            data_dict = {f"{process_name}_{y}{sys}": df}
        return data_dict

    def process_sig(self, data_dict_, y='2024'):
        sig_groups = ['ttZ', 'ttH']
        sig_names = []
        for g in sig_groups:
            #sig_names = sig_names + re.findall(rf'{g}_201\d\w*', ' '.join(data_dict_.keys()))
            sig_names = sig_names + re.findall(rf'{g}_202\d\w*', ' '.join(data_dict_.keys()))
        for sig_name in sig_names:
            s_pre = re.search(r'(ttH|ttZ)',sig_name).group()
            if sig_name not in data_dict_: continue
            #print("made it")
            df = data_dict_[sig_name]
            for i,_ in enumerate(self.pt_bins[:-1]):
                new_sig_name = f"{s_pre}{i}{sig_name.replace(s_pre, '')}"
                data_dict_[new_sig_name] = df[(df['genZHpt']>= self.pt_bins[i]) &
                                              (df['genZHpt'] < self.pt_bins[i+1])]
            data_dict_[f"{s_pre}{len(self.pt_bins) -1}{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            #data_dict_[f"{s_pre}0{sig_name.replace(s_pre, '')}"] = df[df['genZHpt'] >= self.pt_bins[-1]]
            data_dict_.pop(sig_name)
        #print(data_dict_)
            
    #def getZhbbWeight(self, df_, year):
    #    tot_weight = (df_['norm_weight'] * np.sign(df_['genWeight']) * df_['topptWeight'] * 
    #              df_['ele_reco_sf'] * df_['ele_id_sf'] * df_['mu_id_sf'] * df_['mu_iso_sf'] *
    #              df_['bbtag_sf'] * df_['btag_sf'] * df_['puWeight']) # and other weights
    #    return tot_weight

    def initialize_hists(self):
        self.histos = {}
        edges = []
        for s,v in self.data_dict.items():
            print(s)
            #y = re.search(r'201\d',s).group()
            y = re.search(r'202\d',s).group()
            w = getZhbbWeight(v,y) if 'data' not in s else np.ones_like((v[nn] if nn in v else v.iloc [:,0]).to_numpy())
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
        all_mc = self.accepted_sig + self.accepted_bkg
        all_but_ttbb = self.accepted_sig + ['TTBar'] #,'ttX','VJets','single_t']
        tth_sig = [s for s in self.accepted_sig if 'ttH' in s]
        ttz_sig = [s for s in self.accepted_sig if 'ttZ' in s]
        ttbar_mc = ['TTBar','tt_B']
        #jec_mc   = ttbar_mc + ['single_t'] + tth_sig + ttz_sig FIXME
        jec_mc   = ttz_sig
        #
        self.write2dc(f'# Rate uncertainties\n')
        # new lumi
        #Systematic('lumi_13TeV_2016',       'lnN',  all_mc, 1.01)
        #Systematic('lumi_13TeV_2024',       'lnN',  all_mc, 1.02)
        #Systematic('lumi_13TeV_2018',       'lnN',  all_mc, 1.015)
        #Systematic('lumi_13TeV_correlated', 'lnN',  all_mc, {'2016':1.006, '2024': 1.009, '2018': 1.02}) # 0.6, 0.9, 2.0
        #Systematic('lumi_13TeV_1718',       'lnN',  all_mc, {'2024':1.006, '2018': 1.002}) # 0.6, 0.2
        #Systematic('lumi_13TeV_correlated', 'lnN',  all_mc, {'2024': 1.009}) # 0.6, 0.9, 2.0
        #Systematic('lumi_13TeV_1718',       'lnN',  all_mc, {'2024':1.006}) # 0.6, 0.2
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
        # Shape Systematics
        self.write2dc(100*'-'+'\n')
        self.write2dc('# Shape uncertainties \n')
        #ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)#, self.hist3d, self.ptclip)
        # when defining shape, must incluse whether it is a mcsta, scale, or up/down syst
        # correlation is determined via having the same across datacards
        # for y in self.years:
        #     # lepton sfs
        #     self.histos = ShapeSystematic(f'ele_id_sf_{y}', 'shape', 'up/down', all_mc, 1, 'ele_id_sfup','ele_id_sfdown').get_shape()
        #     self.histos = ShapeSystematic(f'ele_reco_sf_{y}', 'shape', 'up/down', all_mc, 1, 'ele_reco_sfup','ele_reco_sfdown').get_shape()
            
        #     # btag sfs
        #     self.histos = ShapeSystematic(f'btag_sf_hf_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sfhf_up','btag_sfhf_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_lf_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sflf_up','btag_sflf_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_hfstats1_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sfhfstats1_up','btag_sfhfstats1_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_hfstats2_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sfhfstats2_up','btag_sfhfstats2_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_lfstats1_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sflfstats1_up','btag_sflfstats1_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_lfstats2_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sflfstats2_up','btag_sflfstats2_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_cferr1_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sfcferr1_up','btag_sfcferr1_down').get_shape()
        #     self.histos = ShapeSystematic(f'btag_sf_cferr2_{y}', 'shape', 'up/down', all_mc, 1, 'btag_sfcferr2_up','btag_sfcferr2_down').get_shape()
        #     self.histos = ShapeSystematic(f'puWeight_{y}',      'shape', 'up/down', all_mc, 1, 'puWeight_up','puWeight_down', extraQC=True).get_shape()
        #     # bbtag sf
        #     self.histos = ShapeSystematic(f'bbtag_sf_{y}', 'shape', 'up/down', all_mc, 1, 'bbtag_sfup','bbtag_sfdown').get_shape()
        #     # jec
        # self.histos = ShapeSystematic(f'jesAbsoluteStat', 'shape', 'qconly', jec_mc, 1, extraQC=True).get_shape()
            
        # #self.histos = ShapeSystematic(f'toppt', 'shape', 'up/down', ttbar_mc, 1, 'topptWeight_Up' ,'topptWeight_Down').get_shape() # using hacky unc. FIXME
        # self.histos = ShapeSystematic(f'mu_id_sf', 'shape', 'up/down', all_mc, 1, 'mu_id_sfup','mu_id_sfdown').get_shape()
        # self.histos = ShapeSystematic(f'mu_iso_sf', 'shape', 'up/down', all_mc, 1, 'mu_iso_sfup','mu_iso_sfdown').get_shape()
        # self.histos = ShapeSystematic(f'mu_reco_sf', 'shape', 'up/down', all_mc, 1, 'mu_reco_sfup','mu_reco_sfdown').get_shape()
        # self.histos = ShapeSystematic(f'isr_tt', 'shape', 'ps', ['TTBar'], 1, 'isr_up','isr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'fsr_tt', 'shape', 'ps', ['TTBar'], 1, 'fsr_up','fsr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'isr_tth', 'shape', 'normup/down', tth_sig, 1, 'isr_up','isr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'isr_ttz', 'shape', 'normup/down', ttz_sig, 1, 'isr_up','isr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'fsr_tth', 'shape', 'normup/down', tth_sig, 1, 'fsr_up','fsr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'fsr_ttz', 'shape', 'normup/down', ttz_sig, 1, 'fsr_up','fsr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic('isr_ttbb', 'shape', 'ps', ['tt_B'], 1, 'isr_up','isr_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic('fsr_ttbb', 'shape', 'ps', ['tt_B'], 1, 'fsr_up','fsr_down', extraQC=True).get_shape()

        # self.histos = ShapeSystematic(f'mu_r_tt', 'shape', 'normup/down', ['TTBar'], 1, 'mu_r_up','mu_r_down').get_shape()
        # self.histos = ShapeSystematic(f'mu_f_tt', 'shape', 'normup/down', ['TTBar'], 1, 'mu_f_up','mu_f_down').get_shape()
        # self.histos = ShapeSystematic(f'mu_r_tth', 'shape', 'normup/down', tth_sig, 1, 'mu_r_up','mu_r_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'mu_f_tth', 'shape', 'normup/down', tth_sig, 1, 'mu_f_up','mu_f_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'mu_r_ttz', 'shape', 'normup/down', ttz_sig, 1, 'mu_r_up','mu_r_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic(f'mu_f_ttz', 'shape', 'normup/down', ttz_sig, 1, 'mu_f_up','mu_f_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic('mu_r_ttbb', 'shape', 'normup/down', ['tt_B'], 1, 'mu_r_up','mu_r_down', extraQC=True).get_shape()
        # self.histos = ShapeSystematic('mu_f_ttbb', 'shape', 'normup/down', ['tt_B'], 1, 'mu_f_up','mu_f_down', extraQC=True).get_shape()
        
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
        #self.write2dc('theo group = CMS_ttbbnorm hdamp_ttbb hdamp UE toppt\n')
        self.write2dc('theo group = CMS_ttbbnorm hdamp_ttbb hdamp UE\n')
        #self.write2dc('theo group += mu_f_ttbb mu_r_ttbb mu_f_tt mu_r_tt mu_f_tth mu_r_tth mu_f_ttz mu_r_ttz\n')
        #self.write2dc('theo group += isr_ttbb fsr_ttbb isr_tt fsr_tt isr_tth fsr_tth isr_ttz fsr_ttz\n')
        #self.write2dc('theo group += tth_ggpdf ttz_ggpdf tth_qsc ttz_qsc ggpdf qqpdf qgpdf tt_qsc ttx_qsc singlet_qsc v_qsc\n')
        #
        
    def setup_Systematics(self):
        process_line = np.array([self.accepted_sig + self.accepted_bkg for _ in range(self.dc_bins)]).flatten()
        Systematic.set_dc_processes(self.dc_dict, process_line)
        ShapeSystematic.set_df_histos_histfuncs(self.data_dict, self.histos)
        
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
            #print(p)
            y = re.findall(r'202\d',p)[0] # first instance of this should be the process year
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
                temp_dict = {'sumw' : to_flat(v['sumw'])}#* (1 if y != '2024' else 3.3032)}#2.2967)} # to just scale to full run2
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
        self.years     = re.findall(r'202\d', name) if name != 'jesHEMIssue' else ['2018']
        if len(self.years) == 0: self.years = ['2024']
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
        cls.p_norms = json.load(open(f'./process_norms/process_norms_run3.json','r'))

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

class ShapeSystematic(Systematic): # Class to handle Datacard shape systematics
    
    '''
    Create additional histograms with systematic variations and adds them to histo dict
    Supports MCStats
    Inheirits from the Systematics Class
    Need to set cut_op and bins_dict to make effective use of this class
    '''

    bins_dict = None
    cut_op   = None
    df       = None

    def __init__(self, name, stype, subtype, process_ids, value, up=None, down=None, info=None, extraQC=False):
        super().__init__(name, stype, process_ids, value, info)
        self.up     = up
        self.down   = down
        self.subtype = subtype
        #self.extraQC = False # to make smoothing before and after plots
        self.extraQC = extraQC
        
    #def getZhbbWeight(self, df_, year):
    #    tot_weight = (df_['norm_weight'] * np.sign(df_['genWeight']) * df_['topptWeight'] * 
    #              df_['ele_reco_sf'] * df_['ele_id_sf'] * df_['mu_id_sf'] * df_['mu_iso_sf'] *
    #              df_['bbtag_sf'] * df_['btag_sf'] * df_['puWeight']) # and other weights
    #    return tot_weight

    def get_shape(self):
        fun_dict = {'mcstat' :self.makeMCStatHist,
                    #'scale'  :self.makeScaleHist,
                    'up/down'    :self.makeUpDownHist,
                    'normup/down':self.makeUpDownHist,
                    'ps'     :self.makeUpDownHist,#self.handlePSUpDown,
                    'pdfas'  :self.makeUpDownPdfasHist,
                    'erdOn'  :self.handleErdOn} # unique case where PS weights for 2024 tth and PS weight for 2016 and 2024 need to be handled 
        if self.subtype in fun_dict: fun_dict[self.subtype]()
        #
        #if False: # testing qc vs no qc 
        if self.extraQC and __name__ == '__main__': self.handle_extraQC()
        self.handleQC()
        #
        return self.histos

    @classmethod
    def set_df_histos_histfuncs(cls, data_dict, hist_dict, 
                                sumw_sumw2=get_sumw_sumw2):# hist_func, ptclip):
        cls.data_dict = data_dict
        cls.histos = hist_dict # to basically add histos dict to scope
        cls.get_sumw_sumw2 = staticmethod(sumw_sumw2)


    def handleErdOn(self):
        for process in self.ids:
            for y in self.years:
                df = self.data_dict[f'{process}_{y}_erdOn']
                weight =  getZhbbWeight(df,y)
                sumw, sumw2 = get_sumw_sumw2(df,w_var,y)
                self.histos[f'{process}_{y}_{self.name}'] = {'sumw':sumw, 
                                                             'sumw2': sumw2}
                

    def handlePSUpDown(self):
        self.makeUpDownHist() # add histos to histos dict
        # special handling for ttHbb[0,1,2,3]_2024 and ttZbb[0,1,2,3]_[2016,2024]
        for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
            for i in range(4):
                sumw_var = ((self.histos[f'ttHbb{i}_2016_{self.name}{ud_str}']['sumw']/self.histos[f'ttHbb{i}_2016']['sumw']) +  # we want the average deviation
                            (self.histos[f'ttHbb{i}_2018_{self.name}{ud_str}']['sumw']/self.histos[f'ttHbb{i}_2018']['sumw']))/2
                sumw2_var = ((self.histos[f'ttHbb{i}_2016_{self.name}{ud_str}']['sumw2']/self.histos[f'ttHbb{i}_2016']['sumw2']) + 
                            (self.histos[f'ttHbb{i}_2018_{self.name}{ud_str}']['sumw2']/self.histos[f'ttHbb{i}_2018']['sumw2']))/2
                self.histos[f'ttHbb{i}_2024_{self.name}{ud_str}'] = {'sumw': np.nan_to_num(sumw_var) * self.histos[f'ttHbb{i}_2024']['sumw'],
                                                                     'sumw2': np.nan_to_num(sumw2_var) * self.histos[f'ttHbb{i}_2024']['sumw2']}

                #sumw_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw']/np.where(self.histos[f'ttZbb{i}_2018']['sumw'] <= 0, 
                #0.00001,self.histos[f'ttZbb{i}_2018']['sumw'])
                sumw_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw']/self.histos[f'ttZbb{i}_2018']['sumw']

                sumw2_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw2']/self.histos[f'ttZbb{i}_2018']['sumw2']
                #sumw2_var = self.histos[f'ttZbb{i}_2018_{self.name}{ud_str}']['sumw2']/np.where(self.histos[f'ttZbb{i}_2018']['sumw2']<=0, 
                #0.00001, self.histos[f'ttZbb{i}_2018']['sumw2'])
                handles_div0 = (lambda x: np.nan_to_num(np.where( ((abs(x)>3.0) | (abs(self.histos[f'ttZbb{i}_2018']['sumw']) <0.05)) ,1, x)))
                self.histos[f'ttZbb{i}_2016_{self.name}{ud_str}'] = {'sumw'  :  handles_div0(sumw_var) * self.histos[f'ttZbb{i}_2016']['sumw'],
                                                                     'sumw2' :  handles_div0(sumw2_var) * self.histos[f'ttZbb{i}_2016']['sumw2']}
                self.histos[f'ttZbb{i}_2024_{self.name}{ud_str}'] = {'sumw'  :  handles_div0(sumw_var) * self.histos[f'ttZbb{i}_2024']['sumw'],
                                                                     'sumw2' :  handles_div0(sumw2_var) * self.histos[f'ttZbb{i}_2024']['sumw2']}
                
    def makeUpDownPdfasHist(self):
        for process in self.ids:
            for y in ['2018','2016','2024']: # special order because of the lack of good shapes for ttH in 2017 #FIXME
                if f'{process}_{y}' not in self.data_dict: continue
                if process == 'ttH' and y == '2024':
                    # take from 2016 and 2018, only for ttH
                    for ud_str in ['Up','Down']:
                        sumw_var = (self.histos[f'{process}_2016_{self.name}{ud_str}']['sumw']/self.histos[f'{process}_2016']['sumw'] + \
                                    self.histos[f'{process}_2018_{self.name}{ud_str}']['sumw']/self.histos[f'{process}_2018']['sumw']) / 2
                        self.histos[f'{process}_{y}_{self.name}{ud_str}'] = {'sumw':sumw_var * self.histos[f'{process}_{y}']['sumw'], 
                                                                             'sumw2':sumw_var * self.histos[f'{process}_{y}']['sumw']}
                #
                df = self.data_dict[f'{process}_{y}']
                nom_weight = getZhbbWeight(df,y)
                pdfas = np.array([self.get_sumw_sumw2(df,nom_weight*df[f'pdfweight_{i}'],y)[0] 
                                 for i in range(103 if process != 'tt_B' else 101)])
                if 'pdf' in self.name:
                    pdfas_type = 'hess' if process != 'tt_B' else 'replica'
                else:
                    pdfas_type = 'alphas'
                pdfas_up, pdfas_down = self.calc_pdfas_unc(pdfas, pdfas_type)
                #
                p_norms_key = re.sub(r'tt_B','ttbb', process) if 'tt_B' in process else process#re.sub(r'\d$','', process)
                
                #sum(weight)/sum(w_var) # to keep the nominal normalization
                pdfas_up   = pdfas_up   * self.p_norms[y][p_norms_key][self.up] 
                pdfas_down = pdfas_down * self.p_norms[y][p_norms_key][self.down] 
                # SO THAT WE ONLY MODEL ACCEPTANCE EFFECTS
                # 
                self.histos[f'{process}_{y}_{self.name}Up'] = {'sumw' :pdfas_up, 
                                                               'sumw2':pdfas_up}
                self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw' :pdfas_down, 
                                                                 'sumw2':pdfas_down}

    @staticmethod
    def calc_pdfas_unc(_pdf, _type='hess'):
        var_pdf = {
            'hess'   : (lambda _x: np.sqrt(np.sum(np.power(_x[1:-2]-_x[0],2), axis=0))),
            'replica': (lambda _x: (np.quantile(_x[1:],.84, axis=0) - np.quantile(_x[1:],.16, axis=0)) / 2),
            'alphas' : (lambda _x: (_x[-2] - _x[-1])/2),
        }
        nom_pdf = {
            'hess'   :(lambda _x: _x[0]),
            'replica':(lambda _x: (np.quantile(_x[1:],.84,axis=0) + np.quantile(_x[1:],.16,axis=0)) / 2),
            'alphas' :(lambda _x: _x[0]),
        }
        err = var_pdf[_type](_pdf)
        nom = nom_pdf[_type](_pdf)
        unc = np.nan_to_num(err/nom, nan=1)
        _out = [_pdf[0]*(1.+unc), _pdf[0]*(1.-unc)]
        return _out # returns up, down
        

        
    def makeUpDownHist(self):
        for process in self.ids:
            for y in self.years:
                if f'{process}_{y}' not in self.data_dict: continue
                df = self.data_dict[f'{process}_{y}']
                #print(df.keys())
                #w_nom  = self.up.replace('_'+self.up.split('_')[-1],'') # format should be weightName_up/down
                w_nom  = self.up.replace(re.search(r'(Up|up)\w*',self.up).group(),'').rstrip('_') # format should be weightName_up/down_foobaropt
                #print(w_nom)
                #print(self.name, process)
                if w_nom == 'pdfWeight' or 'mu_' in w_nom or self.subtype == 'ps' or 'isr' in w_nom or 'fsr' in w_nom: 
                    df[w_nom] = 1.0
                weight =  getZhbbWeight(df,y)
                for ud_str, ud in zip(['Up','Down'],[self.up,self.down]):
                    nominal_weight = weight
                    #if 'new_tt_' in process and self.subtype is 'ps':
                        #nominal_weight = (weight/df['weight']) * df[f'{ud}_weight'] # to insert correct event weight
                    #    if 'mu_' in w_nom: df[w_nom] = df[f'{ud}_weight'] # to fix issue with tt_bb
                    w_var = nominal_weight*df[ud]/df[w_nom]
                    if 'norm' in self.subtype or ('tt_' in process and self.subtype == 'ps'): 
                        #p_norms_key = re.sub(r'tt_','tt', process) if 'tt_bb' in process or 'tt_2b' in process else re.sub(r'\d$','', process)
                        p_norms_key = re.sub(r'tt_B','ttbb', process) if 'tt_B' in process else process#re.sub(r'\d$','', process)\
                        print(self.p_norms[y], ud)
                        w_var = w_var * self.p_norms[y][p_norms_key][ud]#sum(weight)/sum(w_var) # to keep the nominal normalization
                    #
                    #if 'new_tt_' in process and ('mu_' in w_nom or self.subtype is 'ps'):
                    #    print(process, y, ud, ud_str)
                    #    print(sum(weight),sum(nominal_weight),sum(w_var)) # they are equal but u/d nominal are not equal
                    #
                    sumw, sumw2 = self.get_sumw_sumw2(df,w_var,y)
                    self.histos[f'{process}_{y}_{self.name}{ud_str}'] = {'sumw':sumw,
                                                                         'sumw2':sumw2}
                    
                
    def makeMCStatHist(self):
        for process in self.ids:
            for y in self.years:
                nom_hist = self.histos[f'{process}_{y}']
                stat_err = np.sqrt(nom_hist['sumw2'])
                mcstat_up, mcstat_down = nom_hist['sumw']+stat_err, nom_hist['sumw']-stat_err
                self.histos[f'{process}_{y}_{self.name}Up'] =   {'sumw' :mcstat_up,
                                                                 'sumw2':mcstat_up} # dont really care about sumw2 here
                self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw' :mcstat_down,
                                                                 'sumw2':mcstat_down} # dont really care about sumw2 here
                
    def handleQC(self):
        for process in self.ids:
            for y in self.years:
                nom     = np.float32(self.histos[f'{process}_{y}']['sumw'])
                nom_err = np.sqrt(np.float32(self.histos[f'{process}_{y}']['sumw2']))
                #print(f'{process}_{y}_{self.name}Up')
                if re.search(r'(jms)|(jmr)', self.name) is not None and f'{process}_{y}_{self.name}Up' not in self.histos:
                    tmp_name = re.search(r'(jms)|(jmr)', self.name).group()+f'_{y}'
                    self.histos[f'{process}_{y}_{self.name}Up'] = {'sumw':[], 
                                                                   'sumw2':self.histos[f'{process}_{y}_{tmp_name}Up']['sumw2']}     # need to initiate 
                    self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw':[], 
                                                                     'sumw2':self.histos[f'{process}_{y}_{tmp_name}Down']['sumw2']}
                else:
                    tmp_name = self.name
                up   = np.float32(self.histos[f'{process}_{y}_{tmp_name}Up']['sumw'])
                down = np.float32(self.histos[f'{process}_{y}_{tmp_name}Down']['sumw'])
                # Step 1: kill_sys for bins where the stat err is larger than the nominal
                up   = np.where(nom<nom_err, nom, up)   
                down = np.where(nom<nom_err, nom, down)
                # Step 3: one_sided_sys
                #one_sided_sys = (((UpRatio > 1) & (DownRatio > 1)) | ((UpRatio < 1) & (DownRatio < 1)))
                is_onesided = (( (up > nom) & (down > nom) ) | ( (up < nom) & (down < nom) ))
                geo_mean = np.sqrt( up * down )
                up , down = np.where(is_onesided, np.divide(up*nom, geo_mean) , up), np.where(is_onesided, np.divide(down*nom, geo_mean) , down) 
                #
                up, down = np.where((np.isinf(up)) | (np.isnan(up)), nom, up), np.where((np.isinf(down)) | (np.isnan(down)), nom, down) # handle cases where div by zero goes to inf
                
                # save hists 
                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = up     # handle nan later
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = down   # handle nan later in TH1
        
    def handle_extraQC(self):
        for process in self.ids:
            for y in self.years:
                # get variation from merged bins for certain systematic
                # and distribute variation to finer bins, 
                # merged bins should be somewhat correlated in nature
                print(self.histos.keys())
                nom     = np.float32(self.histos[f'{process}_{y}']['sumw'])
                nom_err = np.sqrt(np.float32(self.histos[f'{process}_{y}']['sumw2']))
                #print(f'{process}_{y}_{self.name}Up')
                if re.search(r'(jms)|(jmr)', self.name) is not None and f'{process}_{y}_{self.name}Up' not in self.histos:
                    tmp_name = re.search(r'(jms)|(jmr)', self.name).group()+f'_{y}'
                    self.histos[f'{process}_{y}_{self.name}Up'] = {'sumw':[], 
                                                                   'sumw2':self.histos[f'{process}_{y}_{tmp_name}Up']['sumw2']}     # need to initiate 
                    self.histos[f'{process}_{y}_{self.name}Down'] = {'sumw':[], 
                                                                     'sumw2': self.histos[f'{process}_{y}_{tmp_name}Down']['sumw2']}
                else:
                    tmp_name = self.name
                up   = np.float32(self.histos[f'{process}_{y}_{tmp_name}Up']['sumw'])
                #up   = np.float32(self.histos[f'{process}_{y}_{tmp_name}_Up']['sumw'])
                down = np.float32(self.histos[f'{process}_{y}_{tmp_name}Down']['sumw'])
                #
                # shape (3, 4, 4) pt,nn,sdm

                #if 'JES' in self.name or 'JEC' in self.name: # sum across NN per pt,sdM
                if 'hdamp' not in self.name and 'UE' not in self.name:
                    for i in range(nom.shape[-1]): # for sdM
                        for j in range(nom.shape[0]): # for pt
                            if 'fsr' in self.name : # merge across NN bins
                                up[j,:,i]   = nom[j,:,i] * np.nansum(up[j,:,i])/np.nansum(nom[j,:,i])
                                down[j,:,i] = nom[j,:,i] * np.nansum(down[j,:,i])/np.nansum(nom[j,:,i])
                            else : # dont merge variation in  first 2 bins
                                up[j,2:,i]   = nom[j,2:,i] * np.nansum(up[j,2:,i])/np.nansum(nom[j,2:,i])
                                down[j,2:,i] = nom[j,2:,i] * np.nansum(down[j,2:,i])/np.nansum(nom[j,2:,i])


                else: # dedicated sample sys: hdamp or UE # sum across sdm and NN
                    for j in range(nom.shape[0]): # pt
                        up[j,:,:]   = nom[j,:,:] * np.nansum(up[j,:,:])/np.nansum(nom[j,:,:])
                        down[j,:,:]   = nom[j,:,:] * np.nansum(down[j,:,:])/np.nansum(nom[j,:,:])
                    #up[:,:,:]   = nom[:,:,:] * np.nansum(up[:,:,:])/np.nansum(nom[:,:,:])
                    #down[:,:,:]   = nom[:,:,:] * np.nansum(down[:,:,:])/np.nansum(nom[:,:,:])
                    # ======
                    # dont merge variation in first 2 bins ===== convener question
                    #up[j,2:,:]   = nom[j,2:,:] * np.nansum(up[j,2:,:])/np.nansum(nom[j,2:,:])
                    #down[j,2:,:] = nom[j,2:,:] * np.nansum(down[j,2:,:])/np.nansum(nom[j,2:,:])

                self.histos[f'{process}_{y}_{self.name}Up']['sumw']   = up     # handle nan later
                self.histos[f'{process}_{y}_{self.name}Down']['sumw'] = down   # handle nan later in TH1
        
        
if __name__ == '__main__':
    hist2d = DataCardShapes(pt_bins,mass_bins,isblind=isblind)
    MakeDataCard(isblind=isblind).makeDC()


















    