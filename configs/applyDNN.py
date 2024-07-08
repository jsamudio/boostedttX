import outvars
import pandas as pd
import awkward as ak
import numpy as np
def dnn_cut(df_):
    base_cuts = (
        (df_['n_b_outZH'] >= 2) &
        (df_['ZH_bbvLscore'] >= 0.6) &
        (df_['n_ak4jets']   >= 5)             &
        #( (df_['isEleE']==True) | (df_['isMuonE']==True)) & # pass sim trigger
        #(df_['passNotHadLep'] == 1) & # might add
        (df_['ZH_pt']       >= 200)& # 200
        (df_['MET_pt']      >= 20)            &
        (df_['ZH_M']        >= 50)            &
        (df_['ZH_M']        <= 200)
    )
    return base_cuts
def applyDNN(events, model_file='newgenm_model.h5'):
    from dnn_model import DNN_model

    dnn_vars = outvars.NN_vars
    cutvars = dnn_vars + ['ZH_pt', 'ZH_M', 'MET_pt']
    a = np.array([ak.to_numpy(events[var], allow_missing=True) for var in cutvars]).T
    data = ak.to_numpy(a)
    df_ = pd.DataFrame(data=data, columns = cutvars)
    resetIndex = (lambda df: df.reset_index(drop=True).copy())
    m_info = {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]],
              'other_settings': {'fl_a': [1, 1.5, 1], 'fl_g': 0.25, 'lr_alpha': 0.0002},
              'n_epochs': 150, 'batch_size': 10256}
    dnn = DNN_model(m_info['sequence'],m_info['other_settings'])
    if 'binary' in model_file:
        nn_model = dnn.Build_Binary_Model(len(dnn_vars), load_weights=model_file)#'nn_ttzh_model.h5')
    else:
        nn_model = dnn.Build_Model(len(dnn_vars), load_weights=model_file)#'nn_ttzh_model.h5')          'n_epochs': 150, 'batch_size': 10256}

    # how do I put together the DF?
    base_cuts = dnn_cut(df_)
    #
    pred_df = resetIndex(df_[dnn_vars][base_cuts])
    print(pred_df)
    if len(pred_df) > 0:
        if 'binary' in model_file:
            pred = nn_model.predict(pred_df)
        else:
            pred = nn_model.predict(pred_df)[:,2]
    else:
        pred = 0
    df_nn_name = model_file.replace("model.h5",'')+'NN'
    df_.loc[:,df_nn_name] = -1.
    df_.loc[base_cuts,df_nn_name] = pred
    events['newgenm_NN'] = df_[df_nn_name]
    print(events['newgenm_NN'])
    #print(df_[df_nn_name])
    # only for current main NN
    #if model_file == 'withbbvlnewgenm_model.h5': # old, but may keep
    #if df_nn_name == cfg.nn: # a bit more robust
    if model_file == 'newgenm_model.h5': # new
        explain_model = dnn.Build_New_Model(len(dnn_vars), nn_model) # output size 64
        if len(pred_df) > 0:
            explain_pred = explain_model.predict(pred_df) # shape (n_events, 64)
        else:
            explain_pred = -10* np.ones((1,m_info['sequence'][-2][-1]))
        print(m_info['sequence'][-2][-1])
        for i in range(m_info['sequence'][-2][-1]):
            df_.loc[:,f'NN_{i}'] = -10
            df_.loc[base_cuts,f'NN_{i}'] = explain_pred[:,i]
