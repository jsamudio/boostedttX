from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import zhbbworkflow_spanet
from zhbbworkflow_spanet import ZHbbBaseProcessor
import outvars

from pocket_coffea.parameters.cuts import passthrough

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(zhbbworkflow_spanet)
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
            "jsons": [f"{localdir}/datasets/ttHTobb.json",
                      #f"{localdir}/datasets/ttHToNonbb.json",
                      #f"{localdir}/datasets/TTZToBB.json",
                      f"{localdir}/datasets/TTZToQQ.json",
                      #f"{localdir}/datasets/TTZToLLNuNu.json",
                      #f"{localdir}/datasets/TTbb_Hadronic.json",
                      f"{localdir}/datasets/TTbb_2L2Nu.json",
                      f"{localdir}/datasets/TTbb_SemiLeptonic.json",
                      #f"{localdir}/datasets/TTToHadronic.json",
                      f"{localdir}/datasets/TTTo2L2Nu.json",
                      f"{localdir}/datasets/TTToSemiLeptonic.json"
                ],
            "filter": {
                "samples":  [
                            "ttHTobb",
                            #"ttHToNonbb",
                            "TTZToQQ",
                            #"TTZToLLNuNu",
                            #"TTZToBB",
                            #"TTbb_Hadronic",
                            "TTbb_SemiLeptonic",
                            #"TTbb_2L2Nu",
                            #"TTToHadronic",
                            #"TTTo2L2Nu",
                            "TTToSemiLeptonic",
                ],
                "samples_exclude": [],
                "year": ['2024']
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
                    'TTZToBB' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                    'TTZToQQ' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                    'TTZToLLNuNu' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                    'ttHTobb' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                    'ttHToNonbb' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                }
            },

        workflow = ZHbbBaseProcessor,
        workflow_options = {"parton_jet_min_dR": 0.3,
                            "parton_jet_min_dR_postfsr": 1.0},
                           #"dump_columns_as_arrays_per_chunk": "/cms/data/jsamudio/boosted/boostedttX/configs/SemiLeptonic/output_columns_parton_matching/" },

        skim = [precut],
        preselections = [event_selection],
        #preselections = [fatjet_cut, lep_cut],
        categories = {
            #"jet_cut": [njet_cut],
            #"met_cut":[met_cut],
            "btag_mask": [btag_mask, vetoE, vetoMu, network_cut],
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
            **count_hist(name="nJetGood", coll="JetGood", bins=8, start=0, stop=8),
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
            #"common": {
                "common": {
                    "inclusive": [ColOut("events", ["n_b_inZH"]),
                        ColOut(
                            "Parton",
                            ["pt", "eta", "phi", "mass", "pdgId", "provenance"]
                        ),
                         ColOut(
                            "PartonMatched",
                            ["pt", "eta", "phi","mass", "pdgId", "provenance", "dRMatchedJet"],
                        ),
                        ColOut(
                            "JetGood",
                            ["pt", "eta", "phi", "hadronFlavour", "btagB"],
                        ),
                        # ColOut(
                        #     "JetGoodTruncated",
                        #     ["pt", "eta", "phi", "hadronFlavour", "btagB"],
                        # ),
                        ColOut(
                            "JetGoodMatched",
                            ["pt", "eta", "phi", "hadronFlavour", "btagB", "dRMatchedJet"],
                        ),
                        ColOut("LeptonGood",
                               ["pt","eta","phi", "pdgId", "charge"],
                               pos_end=1, store_size=False),
                        ColOut("PuppiMET", ["phi","pt","significance"]),
                        ColOut("Generator",["x1","x2","id1","id2","xpdf1","xpdf2"]),
                        ColOut("LeptonParton",["pt","eta","phi","mass","pdgId"]),
                        ColOut("FatJetSorted",["pt", "eta", "phi", "mass", "btagBB"], pos_end=1, store_size=False),
                        ColOut(
                            "FatJetMatched",
                            ["pt", "eta", "phi", "mass", "dRMatchedJet", "btagBB"],
                        ),
                        ColOut(
                            "HiggsMatched",
                            ["pt", "provenance"],
                        ),
                ],

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
                #    },
                #"TTZToQQ": {
                #    "bycategory": {
                #        "btag_mask": [
                #            ColOut("HiggsGen",
                #                   ["pt", "eta", "phi", "mass", "pdgId"], pos_end=1, store_size=False),
                #            ]
                #        }
                #    }
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
        "executor"       : "futures",
        "env"            : "myenv",
        "cores"          : 4,
        "workers"        : 1,
        "scaleout"       : 50,
        "queue"          : "microcentury",
        "walltime"       : "00:40:00",
        "mem_per_worker" : "16GB", # GB
        "disk_per_worker" : "16GB", # GB
        "exclusive"      : False,
        "chunk"          : 400000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False
    }
