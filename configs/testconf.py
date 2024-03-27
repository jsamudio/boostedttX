from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import zhbbworkflow
from zhbbworkflow import ZHbbBaseProcessor

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(zhbbworkflow)
cloudpickle.register_pickle_by_value(custom_cut_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Load defaults?

from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/skim_params.yaml",
                                                  update = True)
NN_vars = [
    'matchedGen_ZHbb_bb',
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


cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [f"{localdir}/datasets/ttHTobb_M125.json",
                      f"{localdir}/datasets/TTZToBB.json",
                      f"{localdir}/datasets/TTToSemiLeptonic.json"
                ],
            "filter": {
                "samples": ["ttHTobb", "TTZToBB", "TTToSemiLeptonic"],
                "samples_exclude": [],
                "year": ['2017']
                }
            },

        workflow = ZHbbBaseProcessor,

        skim = [precut],
        preselections = [event_selection],
        categories = {
            "baseline": [],
            "btag_mask": [btag_mask]
            },
        weights = {
            "common": {
                "inclusive": ["genWeight", "lumi", "XS", "pileup"],
            }
        },
        variations = {
            "weights": {
                "common": {
                    "inclusive": ["pileup"],
                }
            }
        },
        variables = {
            #**muon_hists(coll="MuonGood", pos=0),
            #**ele_hists(coll="ElectronGood", pos=0),
            #**count_hist(name="nElectronGood", coll="ElectronGood", bins=3, start=0, stop=3),
            #**count_hist(name="nJetGood", coll="JetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nbJetGood", coll="bJetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nFatJetGood", coll="FatJetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nLeptonGood", coll="LeptonGood", bins=3, start=0, stop=3),
            "mAK8" : HistConf([Axis(coll="FatJetGood", field="particleNet_mass", bins = 100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
            "zhbbtag" : HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            "zhbbtag_sorted" : HistConf([Axis(coll="FatJetSorted", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            "outZH_b1_pt" : HistConf([Axis(coll="events", field="outZH_b1_pt", bins = 100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)])
        },
        columns = {
            "common": {
                "inclusive": [ColOut("events", NN_vars)],
                "bycategory": {}
            },
            "bysample": {
                #"ttHTobb": {"inclusive": [ColOut("events", NN_vars)]},
            }
        }
        )

run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 1,
        "scaleout"       : 50,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
        "queue"          : "microcentury",
        "walltime"       : "00:40:00",
        "mem_per_worker" : "4GB", # GB
        "disk_per_worker" : "1GB", # GB
        "exclusive"      : False,
        "chunk"          : 400000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False
    }
