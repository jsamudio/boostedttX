'''
Lists of output variables
'''

sig_vars = ['process', 'Zbb', 'Hbb', 'Hnonbb', 'Zqq',
        'Zllnunu', 'matchedGenZH', 'matchedGen_Zbb',
        'matchedGen_Hbb', 'matchedGen_ZHbb', 'matchedGen_Zqq',
        'matchedGen_ZHbb_bb', 'matchedGen_ZHbb_b', 'matchedGen_ZHbb_nob',
        'nJetGood', 'ZH_pt', 'MET_pt', 'ZH_M']

bkg_vars = ['process', 'tt_type',
        'nJetGood', 'ZH_pt', 'MET_pt', 'ZH_M']

weight_vars = ['genWeight', 'norm_weight', 'topptWeight']

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
    'n_b_outZH', 'n_q_outZH', "ZH_bbvLscore"]
