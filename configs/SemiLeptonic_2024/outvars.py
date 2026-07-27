'''
Lists of output variables
'''

common_vars = ['process', 'ZH_pt', 'MET_pt', 'ZH_M', 'n_ak4jets', 'n_b_outZH','n_b_inZH', 'n_q_inZH', 'n_q_outZH', 'n_ak4_outZH', 'ZH_bbvLscore', 'ZH_rap',
               'ZH_GloParT3_Xbb', 'ZH_GloParT3_QCD', 'ZH_PNet_XbbVsQCD', 'ZH_PNetLegacy_XbbVsQCD', 
              "encoderFeature1", "encoderFeature2", "encoderFeature3", "encoderFeature4", "encoderFeature5", "encoderFeature6", "encoderFeature7", "encoderFeature8",
              "encoderFeature9", "encoderFeature10", "encoderFeature11", "encoderFeature12", "encoderFeature13", "encoderFeature14", "encoderFeature15", "encoderFeature16",
              "encoderFeature17", "encoderFeature18", "encoderFeature19", "encoderFeature20", "encoderFeature21", "encoderFeature22", "encoderFeature23", "encoderFeature24",
              "encoderFeature25", "encoderFeature26", "encoderFeature27", "encoderFeature28", "encoderFeature29", "encoderFeature30", "encoderFeature31", "encoderFeature32"]

trig_vars = ['Ele30_WPTight_Gsf', 'Ele115_CaloIdVT_GsfTrkIdT', 'Ele50_CaloIdVT_GsfTrkIdT_PFJet165', 'IsoMu24', 'Mu50', 'HighPtTkMu100', 'CascadeMu100']

validation_vars = ['nPV', 'nPVGood', 'MET_pt', 'MET_phi', 'lep_pt', 'ele_pt', 'muon_pt',
                   'lep_eta', 'ele_eta', 'muon_eta', 'n_ak4', 'n_bjet', 'n_ak8',
                   'jet1_pt', 'jet2_pt', 'bjet1_pt', #'bjet2_pt', 
                   'jet1_eta', 'jet2_eta',
                   'bjet1_eta', #'bjet2_eta',
                   'jet1_btag', 'jet2_btag', 'bjet1_btag',
                   #'bjet2_btag',
                   'fatjet1_pt', 'fatjet1_eta', 'fatjet1_mass']
                   
sig_vars = ['Zbb', 'Hbb', 'Hnonbb', 'Zqq',
            'Zllnunu', 'matchedGenZH', 'matchedGen_Zbb',
            'matchedGen_Hbb', 'matchedGen_ZHbb', 'matchedGen_Zqq',
            'matchedGen_ZHbb_bb', 'matchedGen_ZHbb_b', 'matchedGen_ZHbb_nob',
            'genZHpt']
'''
OLD SIG VARS
sig_vars = ['process', 'Zbb', 'Hbb', 'Hnonbb', 'Zqq',
        'Zllnunu', 'matchedGenZH', 'matchedGen_Zbb',
        'matchedGen_Hbb', 'matchedGen_ZHbb', 'matchedGen_Zqq',
        'matchedGen_ZHbb_bb', 'matchedGen_ZHbb_b', 'matchedGen_ZHbb_nob',
        'nJetGood', 'ZH_pt', 'MET_pt', 'ZH_M', 'genZHpt']#, 'ZH_sdm']
'''

bkg_vars = ['tt_B']

wz_bkg_vars = ['nJetGood', 'ZH_pt', 'MET_pt', 'ZH_M'] # probably not needed with common vars now

weight_vars = ['genWeight', 'norm_weight',
               'topptWeight', 'topptWeight_Up', 'topptWeight_Down',
               'ele_reco_sf', 'ele_reco_sfup', 'ele_reco_sfdown',
               'ele_id_sf', 'ele_id_sfup', 'ele_id_sfdown',
               'ele_trig_sf', 'ele_trig_sfup', 'ele_trig_sfdown',
               'mu_id_sf', 'mu_id_sfup', 'mu_id_sfdown', #'mu_id_sfstat',
               'mu_iso_sf', 'mu_iso_sfup', 'mu_iso_sfdown', #'mu_iso_sfstat',
               'mu_trig_sf', 'mu_trig_sfup', 'mu_trig_sfdown', #'mu_iso_sfstat',
               # 'mu_reco_sf', 'mu_reco_sfup', 'mu_reco_sfdown', 'mu_reco_sfstat',
               # #'bbtag_sf', 'bbtag_sfup', 'bbtag_sfdown',
               'btag_sf', 'btag_sfup', 'btag_sfdown',
               # 'btag_sfhf', 'btag_sfhf_up', 'btag_sfhf_down',
               # 'btag_sflf', 'btag_sflf_up', 'btag_sflf_down',
               # 'btag_sfhfstats1', 'btag_sfhfstats1_up', 'btag_sfhfstats1_down',
               # 'btag_sfhfstats2', 'btag_sfhfstats2_up', 'btag_sfhfstats2_down',
               # 'btag_sflfstats1', 'btag_sflfstats1_up', 'btag_sflfstats1_down',
               # 'btag_sflfstats2', 'btag_sflfstats2_up', 'btag_sflfstats2_down',
               # 'btag_sfcferr1', 'btag_sfcferr1_up', 'btag_sfcferr1_down',
               # 'btag_sfcferr2', 'btag_sfcferr2_up', 'btag_sfcferr2_down',
               'puWeight', 'puWeight_up', 'puWeight_down',
               'isr_up', 'isr_down', 'fsr_up', 'fsr_down']
               #'mu_r_up', 'mu_r_down', # scale weight not in some samples
               #'mu_f_up', 'mu_f_down',
               #'mu_rf_up', 'mu_rf_down']

NN_vars = [
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score','outZH_b2_score',
    'outZH_q1_pt','outZH_q2_pt',
    'outZH_q1_score','outZH_q2_score',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1','outZH_q_q_dr_nearb2',
    'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    'outZH_b1_qq_dr','outZH_b2_qq_dr',
    'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr','ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    'l_b2_mtb',
    #
    'ZH_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_ZHbb',
    'outZH_max_ak8pnetMass',
    'outZH_b12_m', 'outZH_b12_dr',
    'ht_b', 'ht_outZH',
    #
    'ak4_bestb_inZH',
    'ak4_worstb_inZH',
    #
    'nonZHbb_q1_dr',
    'nonZHbb_b1_dr',
    'inZHb_outZHb_dr',
    #
    'ZH_l_dr', 'ZH_l_invM',
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZH', 'n_q_inZH',
    'n_b_outZH', 'n_q_outZH', "ZH_bbvLscore"]#, 'min_wpart_ZH_dR', 'max_wpart_ZH_dR', 'n_wpart_ZH_dR_0p4', 'n_wpart_ZH_dR_0p8', 'n_wpart_ZH_dR_1p2', 'min_topb_ZH_dR', 'max_topb_ZH_dR', 'n_topb_ZH_dR_0p4', 'n_topb_ZH_dR_0p8', 'n_topb_ZH_dR_1p2']
    #'ttzbb', 'tthbb', 'ttbb', 'ttlf', 'ttcc']

spanet_vars = [
    'JetGood_pt', 'JetGood_eta', 'JetGood_phi', 'JetGood_btagDeepFlavB', #'JetGood_btagL', 'JetGood_btagM', 'JetGood_btagH', # maybe also mass but they do zeroes_like('pt'), and sin/cos phi
    #'ZH_pt', 'ZH_eta', 'ZH_phi', 'ZH_bbvLscore', # mass decorrelated
    'LeptonGood_pt', 'LeptonGood_eta', 'LeptonGood_phi'] #'LeptonGood_is_electron', # sin and cos phi
    #'MET_pt', 'MET_phi'] #'MET_eta',
    #'ht']
    
