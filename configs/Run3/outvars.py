'''
Lists of output variables of the Run 3 analysis
'''

common_vars = ['process', 'ZH_pt', 'MET_pt', 'ZH_M', 'ZH_rap', 'n_ak4jets',
               'n_b_outZH', 'n_b_inZH', 'n_q_inZH', 'n_q_outZH', 'n_ak4_outZH',
               'ZH_bbvLscore', 'ZH_GloParT3_Xbb', 'ZH_GloParT3_QCD', 'ZH_PNet_XbbVsQCD', 'ZH_PNetLegacy_XbbVsQCD',
               'xbb_genmatch', 'n_genb_matchZH']

# common_vars that only exist for MC (dropped from the data columns)
mc_only_vars = ['xbb_genmatch', 'n_genb_matchZH']

validation_vars = ['nPV', 'nPVGood', 'MET_pt', 'MET_phi', 'lep_pt', 'ele_pt', 'muon_pt',
                   'lep_eta', 'ele_eta', 'muon_eta', 'n_ak4', 'n_bjet', 'n_ak8',
                   'jet1_pt', 'jet2_pt', 'bjet1_pt',
                   'jet1_eta', 'jet2_eta', 'bjet1_eta',
                   'jet1_btag', 'jet2_btag', 'bjet1_btag',
                   'fatjet1_pt', 'fatjet1_eta', 'fatjet1_mass']

sig_vars = ['Zbb', 'Hbb', 'Hnonbb', 'Zqq',
            'Zllnunu', 'matchedGenZH', 'matchedGen_Zbb',
            'matchedGen_Hbb', 'matchedGen_ZHbb', 'matchedGen_Zqq',
            'matchedGen_ZHbb_bb', 'matchedGen_ZHbb_b', 'matchedGen_ZHbb_nob',
            'genZHpt']

bkg_vars = ['tt_B']

# Per-event weights (MC). Trigger SFs and b-tag SFs are 1 until they are derived (see README).
weight_vars = ['genWeight', 'norm_weight',
               'topptWeight', 'topptWeight_Up', 'topptWeight_Down',
               'ele_reco_sf', 'ele_reco_sfup', 'ele_reco_sfdown',
               'ele_id_sf', 'ele_id_sfup', 'ele_id_sfdown',
               'ele_trig_sf', 'ele_trig_sfup', 'ele_trig_sfdown',
               'mu_id_sf', 'mu_id_sfup', 'mu_id_sfdown',
               'mu_iso_sf', 'mu_iso_sfup', 'mu_iso_sfdown',
               'mu_trig_sf', 'mu_trig_sfup', 'mu_trig_sfdown',
               'btag_sf', 'btag_sfup', 'btag_sfdown',
               'puWeight', 'puWeight_up', 'puWeight_down',
               'isr_up', 'isr_down', 'fsr_up', 'fsr_down',
               'dctr_hdamp_up', 'dctr_hdamp_down',
               'mu_r_up', 'mu_r_down', 'mu_f_up', 'mu_f_down', 'mu_rf_up', 'mu_rf_down']

# Per-event factors that enter the nominal event weight on top of norm_weight * sign(genWeight)
nominal_weight_factors = ['topptWeight', 'ele_reco_sf', 'ele_id_sf', 'ele_trig_sf',
                          'mu_id_sf', 'mu_iso_sf', 'mu_trig_sf', 'btag_sf', 'puWeight']

# Trigger-efficiency columns (TrigEff.ipynb): lepton kinematics, event variables, weights; the
# HLT bits are exported from the HLT collection (custom_cut_functions.hlt_paths)
trigeff_vars = ['ele_pt', 'ele_eta', 'muon_pt', 'muon_eta', 'n_ak4jets', 'n_b_outZH']
