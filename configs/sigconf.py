from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import testingworkflow
from testingworkflow import ZHbbTestProcessor
import outvars
from pocket_coffea.parameters.cuts import passthrough

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(testingworkflow)
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
                                                  f"{localdir}/params/sample_params.yaml",
                                                  f"{localdir}/params/event_flags.yaml",
                                                  update = True)

cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [#f"{localdir}/datasets/ttHTobb_M125.json",
                      #f"{localdir}/datasets/ttHToNonbb_M125.json",
                      f"{localdir}/datasets/TTZToBB.json",
                      #f"{localdir}/datasets/TTZToQQ.json",
                      #f"{localdir}/datasets/TTZToLLNuNu.json",
                      #f"{localdir}/datasets/TTbb_Hadronic.json",
                      #f"{localdir}/datasets/TTbb_2L2Nu.json",
                      #f"{localdir}/datasets/TTbb_SemiLeptonic.json",
                      #f"{localdir}/datasets/TTToHadronic.json",
                      #f"{localdir}/datasets/TTTo2L2Nu.json",
                      #f"{localdir}/datasets/TTToSemiLeptonic.json"
                ],
            "filter": {
                "samples":  [
                            #"ttHTobb",
                            #"ttHToNonbb",
                            #"TTZToQQ",
                            #"TTZToLLNuNu",
                            "TTZToBB",
                            #"TTbb_Hadronic",
                            #"TTbb_SemiLeptonic",
                            #"TTbb_2L2Nu",
                            #"TTToHadronic",
                            #"TTTo2L2Nu",
                            #"TTToSemiLeptonic",
                ],
                "samples_exclude": [],
                "year": ['2017']
                },
                "subsamples": {
                }
            },

        workflow = ZHbbTestProcessor,

        skim = [precut],
        preselections = [event_selection],
        categories = {
            "btag_mask": [btag_mask],
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
        save_skimmed_files = "/cms/data/store/user/jsamudio/NanoAODv9/pocketCoffea/",
        variables = {
            #**muon_hists(coll="MuonGood", pos=0),
            #**ele_hists(coll="ElectronGood", pos=0),
            #**count_hist(name="nElectronGood", coll="ElectronGood", bins=3, start=0, stop=3),
            #**count_hist(name="nJetGood", coll="JetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nbJetGood", coll="bJetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nFatJetGood", coll="FatJetGood", bins=8, start=0, stop=8),
            #**count_hist(name="nLeptonGood", coll="LeptonGood", bins=3, start=0, stop=3),
            #"mAK8" : HistConf([Axis(coll="FatJetGood", field="particleNet_mass", bins = 100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
            #"zhbbtag" : HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            #"zhbbtag_sorted" : HistConf([Axis(coll="FatJetSorted", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            #"newgenm_NN" : HistConf([Axis(coll="events", field="newgenm_NN", bins = [0., 0.08306063, 0.43137971, 0.55986929, 0.73463416, 0.8649936, 1. ], start=0, stop=1, label=r"$DNN Score$", pos=0, underflow=False, overflow=False)]),
            #"outZH_b1_pt" : HistConf([Axis(coll="events", field="outZH_b1_pt", bins = 100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)])
        },
        columns = {
            "common": {
                "inclusive": [ColOut("events", ["run", "event", "nFatJetGood", "nElectronGood", "nJetGood", "nMuonGood"])],
                "bycategory": {}
            },
            "bysample": {
                #"ttHTobb": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars)]},
                #"ttHToNonbb": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars)]},
                #"TTZToQQ": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars)]},
                #"TTZToLLNuNu": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars)]},
                #"TTZToBB": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars)]},
                #"TTbb_Hadronic": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
                #"TTbb_SemiLeptonic": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
                #"TTbb_2L2Nu": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
                #"TTToHadronic": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
                #"TTTo2L2Nu": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
                #"TTToSemiLeptonic": {"inclusive": [ColOut("events", ['tt_B']+outvars.NN_vars+outvars.bkg_vars)]},
            }
        }
        )

run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "myenv",
        "cores"          : 4,
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
