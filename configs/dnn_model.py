## temp file to store dnn model
import pandas as pd
import numpy as np
import os
import sys
#if __name__ == '__main__':
#    import subprocess as sb
#    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as K

#import config.ana_cff as cfg

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

nodak8md_dnn_ZH_vars = [
    # may have to get rid of due to bad p-score:
    # n_q_outZh, l_b2_dr, n_ak4jets
    'outZH_b1_pt','outZH_b2_pt',
    'outZH_b1_score','outZH_b2_score',
    'outZh_q1_pt','outZh_q2_pt',
    'outZh_q1_btag','outZh_q2_btag',
    #
    'outZH_b1_q_mindr','outZH_b2_q_mindr',
    'outZH_q_q_dr_nearb1','outZH_q_q_dr_nearb2',
    'outZH_qq_M_nearb1','outZH_qq_M_nearb2',
    #'outZH_b1q_M', bad p-score!
    #'outZH_b2q_M',
    'outZH_b1_qq_dr','outZH_b2_qq_dr',
    'outZH_b1qq_M','outZH_b2qq_M',
    'ZH_b1qq_dr','ZH_b2qq_dr',
    'ZH_lbb1qq_dr','ZH_lbb2qq_dr',
    'l_b2_mtb',
    #
    'Zh_closeb_invM',#'Zh_closeq_invM',
    'n_ak8jets', 'n_ak4jets','n_ak8_Zhbb',
    'outZh_max_ak8sdM',
    'outZh_b12_m', 'outZh_b12_dr',
    'ht_b', 'ht_outZh',
    #
    'ak4_bestb_inZH',
    'ak4_worstb_inZH',
    #
    'nonZhbb_q1_dr',
    'nonZhbb_b1_dr',
    'inZhb_outZhb_dr',
    #
    'Zh_l_dr', 'Zh_l_invM_sd',
    'l_b1_invM','l_b2_invM',
    'l_b1_dr','l_b2_dr',
    #
    'spher','aplan',
    'n_b_inZh', 'n_q_inZh',
    'n_b_outZh', 'n_q_outZh'
]


withbbvl_dnn_ZHgenm_vars = nodak8md_dnn_ZH_vars+['Zh_bbvLscore']


class DNN_model:

    #alpha  = cfg.dnn_ZH_alpha
    #outDir = cfg.DNNoutputDir
    #outName= cfg.DNNoutputName
    #useW   = cfg.DNNuseWeights
    #

    seq_dict = {
        'Dense'  : (lambda x: keras.layers.Dense(x, activation='relu', )),#kernel_regularizer=y(z))),
        'Dropout': (lambda x: keras.layers.Dropout(x)),
        'Batch'  : (lambda : keras.layers.BatchNormalization())
    }

    #@classmethod
    #def focal_loss(cls,y_true, y_pred):
    #    gamma = cls.fl_g
    #    alpha = cls.fl_a
    #    pt_1 = tf.where(tf.equal(y_true, 1), y_pred, tf.ones_like(y_pred))
    #    pt_0 = tf.where(tf.equal(y_true, 0), y_pred, tf.zeros_like(y_pred))
    #    return -K.sum(alpha * K.pow(1. - pt_1, gamma) * K.log(pt_1))-K.sum((1-alpha) * K.pow( pt_0, gamma) * K.log(1. - pt_0))

    def __init__(self, sequence, other_settings):
        self.sequence = sequence # tuple of ('Dense', (128,keras.regularizers.l1,0.0))
        self.fl_alpha = other_settings['fl_a']
        self.fl_gamma = other_settings['fl_g']
        self.lr_alpha = other_settings['lr_alpha']


    @staticmethod
    def cfl(alpha, gamma=2.):
        alpha = np.array(alpha, dtype=np.float32)

        def categorical_focal_loss_fixed(y_true, y_pred):
            """
            :param y_true: A tensor of the same shape as `y_pred`
            :param y_pred: A tensor resulting from a softmax
            :return: Output tensor.
            """
            # Clip the prediction value to prevent NaN's and Inf's
            epsilon = K.epsilon()
            y_pred = K.clip(y_pred, epsilon, 1. - epsilon)
            # Calculate Cross Entropy
            print(y_pred.shape, y_true.shape)
            cross_entropy = -y_true * K.log(y_pred)
            # Calculate Focal Loss
            loss = alpha * K.pow(1 - y_pred, gamma) * cross_entropy
            # Compute mean loss in mini_batch
            return K.mean(K.sum(loss, axis=-1))
            #
        return categorical_focal_loss_fixed


    @staticmethod
    def resetIndex(df_):
        return df_.reset_index(drop=True).copy()


    def Build_Model(self, input_shape, load_weights=None):
        #
        main_input = keras.layers.Input(shape=[input_shape], name='input')
        layer = keras.layers.BatchNormalization()(main_input)
        for seq in self.sequence:
            layer = self.seq_dict[seq[0]](seq[1])(layer)
        #
        output     = keras.layers.Dense(3, activation='softmax', name='output')(layer)
            #
        model      = keras.models.Model(inputs=main_input, outputs=output, name='model')
        optimizer  = keras.optimizers.Adam(learning_rate=self.lr_alpha, clipnorm=1)

        if (load_weights and os.path.exists('./'+load_weights)):
            model.load_weights('./'+load_weights)
        #
        #from focal_loss import BinaryFocalLoss
        model.compile(
            loss=[self.cfl(self.fl_alpha,self.fl_gamma)],
            optimizer=optimizer,
            #optimizer='adam',
            #
            metrics=['accuracy',tf.keras.metrics.AUC()])
        ##
        return model

    def Build_Binary_Model(self, input_shape, load_weights=None):
        #
        main_input = keras.layers.Input(shape=[input_shape], name='input')
        layer = keras.layers.BatchNormalization()(main_input)
        for seq in self.sequence:
            layer = self.seq_dict[seq[0]](seq[1])(layer)
        #
        output     = keras.layers.Dense(1, activation='sigmoid', name='output')(layer)
            #
        model      = keras.models.Model(inputs=main_input, outputs=output, name='model')
        optimizer  = keras.optimizers.Adam(learning_rate=self.lr_alpha)

        if (load_weights and os.path.exists('./'+load_weights)):
            model.load_weights('./'+load_weights)
        #
        #from focal_loss import BinaryFocalLoss
        model.compile(
            loss='binary_crossentropy',
            optimizer=optimizer,
            #optimizer='adam',
            #
            metrics=['accuracy',tf.keras.metrics.AUC()])
        ##
        return model

    def Build_New_Model(self, input_shape, model):
        main_input = keras.layers.Input(shape=[input_shape], name='input')
        layer = keras.layers.BatchNormalization()(main_input)
        for seq in self.sequence[:-2]:
            layer = self.seq_dict[seq[0]](seq[1])(layer)
        #
        output     = keras.layers.Dense(64, name='output')(layer)
        #
        new_model      = keras.models.Model(inputs=main_input, outputs=output, name='new_model')
        optimizer  = keras.optimizers.Adam(learning_rate=self.lr_alpha)

        new_model.set_weights(model.get_weights()[:-2])
        #
        #from focal_loss import BinaryFocalLoss
        new_model.compile(
            loss=[self.cfl(self.fl_alpha,self.fl_gamma)],
            optimizer=optimizer,
            #optimizer='adam',
            #
            metrics=['accuracy',tf.keras.metrics.AUC()])
        ##
        return new_model


def train_model(m_info):
    trainX, trainY, valX, valY, testX, testY, model = prep_model_data(m_info)
    # --
    # earlystopping here?
    cb = [tf.keras.callbacks.EarlyStopping(monitor='val_loss',patience=3,restore_best_weights=True)]
    # --
    print(model.summary())
    # --
    history = model.fit(
        trainX,
        trainY,
        epochs     = m_info['n_epochs'],
        batch_size = m_info['batch_size'],
        validation_data = (valX, valY),
        callbacks = cb,
        shuffle = True,
        verbose = 1,
    )
    if __name__ == '__main__':
        loss ,acc, auc             = model.evaluate(testX, testY)#, sample_weight = testW['DNNweight'].values)
        tr_loss ,tr_acc, tr_auc    = model.evaluate(trainX,trainY)
        val_loss ,val_acc, val_auc = model.evaluate(valX,  valY)
        print("\n\n")
        print(f"Train Set Loss: {tr_loss:10.4f}  Train Set Acc: {tr_acc:10.4f} Train Set AUC: {tr_auc:10.4f}")
        print(f"Val   Set Loss: {val_loss:10.4f}  Val   Set Acc: {val_acc:10.4f} Val   Set AUC: {val_auc:10.4f}")
        print(f"Test  Set Loss: {loss:10.4f}  Test  Set Acc: {acc:10.4f} Test  Set AUC: {auc:10.4f}\n\n")
        plot_history(history)
        print('here')
        # dont overwrite ~~~ !!!
        #exit()
        model.save_weights('./'+'/'+'newgenm_model'+'.weights.h5')
        # ======= testing binary
        #model.save_weights(cfg.dnn_ZH_dir+'/'+'binary_model'+'.h5')
        #model.save_weights(cfg.dnn_ZH_dir+'/'+'newreduced1p0_model'+'.h5')
        print('here')

    return model, testX, testY

def train_binary_model(m_info):
    trainX, trainY, valX, valY, testX, testY, model = prep_model_data(m_info, is_binary=True)
    # --
    # earlystopping here?
    cb = [tf.keras.callbacks.EarlyStopping(monitor='val_loss',patience=3,restore_best_weights=True)]
    # --
    print(model.summary())
    # --
    history = model.fit(
        trainX,
        trainY[:,2],
        epochs     = 40,
        batch_size = m_info['batch_size'],
        validation_data = (valX, valY[:,2]),
        callbacks = cb,
        shuffle = True,
        verbose = 1,
    )
    if __name__ == '__main__':
        loss ,acc, auc             = model.evaluate(testX, testY[:,2])#, sample_weight = testW['DNNweight'].values)
        tr_loss ,tr_acc, tr_auc    = model.evaluate(trainX,trainY[:,2])
        val_loss ,val_acc, val_auc = model.evaluate(valX,  valY[:,2])
        print("\n\n")
        print(f"Train Set Loss: {tr_loss:10.4f}  Train Set Acc: {tr_acc:10.4f} Train Set AUC: {tr_auc:10.4f}")
        print(f"Val   Set Loss: {val_loss:10.4f}  Val   Set Acc: {val_acc:10.4f} Val   Set AUC: {val_auc:10.4f}")
        print(f"Test  Set Loss: {loss:10.4f}  Test  Set Acc: {acc:10.4f} Test  Set AUC: {auc:10.4f}\n\n")
        plot_history(history)
        print('here')
        model.save_weights('./'+'/'+'binary_model'+'.h5')
        #model.save_weights(cfg.dnn_ZH_dir+'/'+'newreduced1p0_model'+'.h5')
        print('here')

    return model, testX, testY

def plot_history(history):
    import matplotlib.pyplot as plt
    hist = pd.DataFrame(history.history)
    hist['epoch'] = history.epoch
    fig, [loss, acc] = plt.subplots(1, 2, figsize=(12, 6))
    loss.set_xlabel('Epoch')
    loss.set_ylabel('Loss')
    loss.grid(True)
    loss.plot(hist['epoch'], hist['loss'],
              label='Train Loss')
    loss.plot(hist['epoch'], hist['val_loss'],
              label = 'Val Loss')
    #loss.set_yscale('log')
    loss.legend

    acc.set_xlabel('Epoch')
    acc.set_ylabel('Acc')
    acc.grid(True)
    acc.plot(hist['epoch'], hist['accuracy'],
                     label='Train Acc')
    acc.plot(hist['epoch'], hist['val_accuracy'],
                     label = 'Val Acc')
    #acc.set_yscale('log')
    acc.legend()
    plt.show()
    plt.close()

def local_test(m_info, train_binary=False):
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import roc_auc_score
    import matplotlib.pyplot as plt
    if train_binary == False:
        model, testX, testY = train_model(m_info)
    else:
        model, testX, testY = train_binary_model(m_info)
    y_pred = model.predict(testX)
    weight = np.ones_like(y_pred[:,2])*.001
    weight = np.where(testY[:,0]==1,.10,weight)
    weight = np.where(testY[:,1]==1,.15,weight)
    print('AUC score', roc_auc_score(testY[:,2],y_pred[:,2],sample_weight=weight))
    cm = confusion_matrix(np.argmax(testY,axis=1), np.argmax(y_pred,axis=1))
    print ("\nThe confusion matrix of the test set on the trained nerual network:\n" , cm)
    # contruct score
    zh_score = cm[:,2]
    print(cm[:,2])
    s_b_val = zh_score[2]*.001/((zh_score[0]*0.10+zh_score[1]*.15)**(1/2))
    s_ttbb_val = zh_score[2]*.001/((zh_score[1]*.15)**(1/2))
    cm = confusion_matrix(np.argmax(testY[y_pred[:,2]>0.6],axis=1), np.argmax(y_pred[y_pred[:,2]>0.6],axis=1))
    print ("\nThe confusion matrix of the test set (high confidence in zh):\n" , cm)
    #out_name = f'{s_b_val:.2f}_{s_ttbb_val:.2f}_{args.job_number}'
    #print(out_name)
    #model.save_weights(cfg.dnn_ZH_dir+'/test_archs/'+out_name+'.h5')

    plt.hist(y_pred[testY[:,2] == 1][:,2],
             bins=10, range=(0,1), histtype='step',
             weights=np.ones(len(y_pred[testY[:,2] == 1][:,2]))*.001, label='SIG')
    plt.hist(y_pred[testY[:,0] == 1][:,2],
             bins=10, range=(0,1), histtype='step',
             weights=np.ones(len(y_pred[testY[:,0] == 1][:,2]))*.10, label='tt')
    plt.hist(y_pred[testY[:,1] == 1][:,2],
             bins=10, range=(0,1), histtype='step',
             weights=np.ones(len(y_pred[testY[:,1] == 1][:,2]))*.15, label='ttbb')
    plt.legend()
    plt.yscale('log')
    plt.xlim(0,1)
    #plt.title(out_name)
    plt.show()
    plt.savefig('nn_training.pdf')


def prep_model_data(m_info, is_binary=False):
    resetIndex = (lambda df: df.reset_index(drop=True).copy())
    #
    trainXY = pd.read_pickle('./trainXY.pkl')
    #trainXY = pd.read_pickle('/cms/data/jsamudio/ttZh_analysis/ttZ-h_eft/DeepSleep/data/NN_files/trainXY.pkl')
    #assert not np.any(np.isnan(trainXY))
    print(trainXY.keys())
    print(trainXY.isna().sum().head(50))
    # get val from trainXY
    valXY   = trainXY.sample(frac=.25, random_state=1)
    trainXY = trainXY.drop(valXY.index).copy()
    #
    testXY  = pd.read_pickle('./testXY.pkl')
    #testXY  = pd.read_pickle('/cms/data/jsamudio/ttZh_analysis/ttZ-h_eft/DeepSleep/data/NN_files/testXY.pkl')
    trainY, valY, testY = [resetIndex(df['label']) for df in [trainXY, valXY, testXY]]
    print(f"Size of : {'train':10} {'val':10} {'test':10}")
    print(f"          {len(trainY):10} {len(valY):10} {len(testY):10}")
    trainX, valX, testX = [resetIndex(df.drop(columns=['label'])) for df in [trainXY, valXY, testXY]]
    print(trainX.shape)
    #
    m_class = DNN_model(m_info['sequence'],m_info['other_settings'])
    #model = m_class.Build_Model(len(cfg.withbbvl_dnn_ZHgenm_vars), )#load_weights='ttzh_model.h5')
    dnn_vars = NN_vars
    #dnn_vars = withbbvl_dnn_ZHgenm_vars
    #dnn_vars = cfg.reduced1p0genm_vars
    #dnn_vars = cfg.withbbvl_dnn_ZH_1p5vars
    # ===== testing binary model
    if is_binary:
        model = m_class.Build_Binary_Model(len(dnn_vars), )#load_weights='ttzh_model.h5')
    else:
        model = m_class.Build_Model(len(dnn_vars), )#load_weights='ttzh_model.h5')
    return (
        trainX.filter(items=dnn_vars).to_numpy(), np.stack(trainY.values), valX.filter(items=dnn_vars).to_numpy(), np.stack(valY.values), testX.filter(items=dnn_vars).to_numpy(), np.stack(testY.values), model
    )


if __name__ == "__main__":
    import json
    import sys
    json_dir = f'./log/nn/'

    # === old NN
    #m_info = {"sequence": [["Dense", 128], ["Dense", 64], ["Dropout", 0.5]], "other_settings": {"fl_a": [0.75, 1, 0.25], "fl_g": 0.5, "lr_alpha": 0.0003}, "n_epochs": 140, "batch_size": 10256}
    # === new NN
    #m_info = {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 'other_settings': {'fl_a': [1, 2, 0.75], 'fl_g': 0.25, 'lr_alpha': 0.0003}, 'n_epochs': 150, 'batch_size': 10256}
    m_info = {'sequence': [['Dense', 128], ['Dense', 64], ['Dropout', 0.5]], 'other_settings': {'fl_a': [2, 2.5, 1], 'fl_g': 0.25, 'lr_alpha': 0.0002}, 'n_epochs': 250, 'batch_size': 10256}

    local_test(m_info)


