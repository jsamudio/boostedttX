from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import spanet_inference_noCuts
from spanet_inference_noCuts import ZHbbSpanetProcessor
import outvars
from pocket_coffea.parameters.cuts import passthrough

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(spanet_inference_noCuts)
cloudpickle.register_pickle_by_value(custom_cut_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

spanet_model_path = "/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalanced24.onnx"

# Load defaults?
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/lepton_scale_factors.yaml",
                                                  f"{localdir}/params/skim_params.yaml",
                                                  f"{localdir}/params/sample_params.yaml",
                                                  f"{localdir}/params/event_flags.yaml",
                                                  update = True)

cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [f"{localdir}/datasets/ttHTobb_M125.json",
                      f"{localdir}/datasets/ttHToNonbb_M125.json",
                      f"{localdir}/datasets/TTZToBB.json",
                      f"{localdir}/datasets/TTZToQQ.json",
                      f"{localdir}/datasets/TTZToLLNuNu.json",
                      f"{localdir}/datasets/TTbb_Hadronic.json",
                      f"{localdir}/datasets/TTbb_2L2Nu.json",
                      f"{localdir}/datasets/TTbb_SemiLeptonic.json",
                      f"{localdir}/datasets/TTToHadronic.json",
                      f"{localdir}/datasets/TTTo2L2Nu.json",
                      f"{localdir}/datasets/TTToSemiLeptonic.json"
                ],
            "filter": {
                "samples":  [
                            #"ttHTobb",
                            #"ttHToNonbb",
                            #"TTZToQQ",
                            #"TTZToLLNuNu",
                            #"TTZToBB",
                            "TTbb_Hadronic",
                            "TTbb_SemiLeptonic",
                            "TTbb_2L2Nu",
                            #"TTToHadronic",
                            #"TTTo2L2Nu",
                            #"TTToSemiLeptonic",
                ],
                "samples_exclude": [],
                "year": ['2017']
                },
                "subsamples": {
                    'TTbb_SemiLeptonic': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                    'TTbb_Hadronic': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                    'TTbb_2L2Nu': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                    'TTToSemiLeptonic': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                    'TTToHadronic': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                    'TTTo2L2Nu': {
                        'tt+LF'   : [get_genTtbarId_100_eq(0)],
                        'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                        'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    },
                }
            },

        workflow = ZHbbSpanetProcessor,
        workflow_options = {"parton_jet_min_dR": 0.3,
                            "parton_jet_min_dR_postfsr": 1.0,
                            "spanet_model" : spanet_model_path},
                           #"dump_columns_as_arrays_per_chunk": "/cms/data/jsamudio/boosted/boostedttX/configs/SemiLeptonic/output_columns_parton_matching/" },

        skim = [passthrough],
        preselections = [passthrough],
        categories = {
            "btag_mask": [passthrough],
            },
        weights = {
            "common": {
                "inclusive": [],
            }
        },
        variations = {
            "weights": {
                "common": {
                    "inclusive": [],
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
            "LHEHT" : HistConf([Axis(coll="LHE", field="HT", bins = 100, start=0, stop=5000, label=r"$HT_{LHE}$ [GeV]")]),
            #"fullHT" : HistConf([Axis(coll="events", field="ttbb_full_ht", bins = 100, start=0, stop=2000, label=r"$HT_{ttbb}$ [GeV]")]),
            "bbHT" : HistConf([Axis(coll="events", field="bb_ht", bins = 100, start=0, stop=4000, label=r"$HT_{bb}$ [GeV]")]),
            "ttHT" : HistConf([Axis(coll="events", field="tt_ht", bins = 100, start=0, stop=4000, label=r"$HT_{tt}$ [GeV]")]),
            "ttbbHT" : HistConf([Axis(coll="events", field="ttbb_full_ht", bins = 100, start=0, stop=5000, label=r"$HT_{ttbb}$ [GeV]")]),
            "ttbbHT_500_750" : HistConf([Axis(coll="events", field="ttbb_full_ht", bins = 100, start=500, stop=750, label=r"$HT_{ttbb}$ [GeV]", underflow=False, overflow=False)]),
            "ttbbHT_750_5000" : HistConf([Axis(coll="events", field="ttbb_full_ht", bins = 100, start=750, stop=5000, label=r"$HT_{ttbb}$ [GeV]", underflow=False, overflow=False)]),
            #"zhbbtag_sorted" : HistConf([Axis(coll="FatJetSorted", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            #"newgenm_NN" : HistConf([Axis(coll="events", field="newgenm_NN", bins = [0., 0.08306063, 0.43137971, 0.55986929, 0.73463416, 0.8649936, 1. ], start=0, stop=1, label=r"$DNN Score$", pos=0, underflow=False, overflow=False)]),
            #"spanet_tthbb" : HistConf(
            #    [Axis(coll="spanet_output", field="tthbb", bins=50, start=0, stop=1, label="tthbb SPANet score")],
            #),
            #"spanet_ttbb" : HistConf(
            #    [Axis(coll="spanet_output", field="ttbb", bins=50, start=0, stop=1, label="ttbb SPANet score")],
            #),
            #"spanet_ttcc" : HistConf(
            #    [Axis(coll="spanet_output", field="ttcc", bins=50, start=0, stop=1, label="ttcc SPANet score")],
            #),
            #"spanet_ttlf" : HistConf(
            #    [Axis(coll="spanet_output", field="ttlf", bins=50, start=0, stop=1, label="ttlf SPANet score")],
            #),
            #"outZH_b1_pt" : HistConf([Axis(coll="events", field="outZH_b1_pt", bins = 100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)])
        },
        columns = {
            "common": {
                "inclusive": [],
                
                "bycategory": {}
            },
            "bysample": {
                #"ttHTobb": {
                #    "bycategory": {
                #        "btag_mask": [
                #            ColOut("HiggsGen",
                #                   ["pt", "eta", "phi", "mass", "pdgId"], pos_end=1, store_size=False),
                #            ]
                #        }
                #,
                #"TTZToBB": {
                #    "bycategory": {
                #        "btag_mask": [
                #            ColOut("HiggsGen",
                #                   ["pt", "eta", "phi", "mass", "pdgId"], pos_end=1, store_size=False),
                #            ]
                #        }
                #    },
                #"ttHTobb": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars+["genZHpt"])]},
                #"ttHToNonbb": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars+["genZHpt"])]},
                #"TTZToQQ": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars+["genZHpt"])]},
                #"TTZToLLNuNu": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars+["genZHpt"])]},
                #"TTZToBB": {"inclusive": [ColOut("events", outvars.NN_vars+outvars.sig_vars+["genZHpt"])]},
                #"TTbb_Hadronic": {"inclusive": [ColOut("events", ["ttbb_full_ht", "ttbb_had_ht", "tt_ht", "bb_ht", "LHE_HT"])]},
                #"TTbb_SemiLeptonic": {"inclusive": [ColOut("events", ["ttbb_full_ht", "ttbb_had_ht", "tt_ht", "bb_ht", "LHE_HT"])]},
                #"TTbb_2L2Nu": {"inclusive": [ColOut("events", ["ttbb_full_ht", "ttbb_had_ht", "tt_ht", "bb_ht", "LHE_HT"])]},
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
