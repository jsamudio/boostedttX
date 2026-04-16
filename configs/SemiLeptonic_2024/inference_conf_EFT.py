from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import spanet_inference
from spanet_inference import ZHbbSpanetProcessor
from pocket_coffea.lib.weights.common.common import common_weights
import outvars
from pocket_coffea.parameters.cuts import passthrough

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(spanet_inference)
cloudpickle.register_pickle_by_value(custom_cut_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))
print(localdir)

from pocket_coffea.lib.calibrators.common.common import default_calibrators_sequence

spanet_model_path = "/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalanced24.onnx"

# Load defaults?
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/lepton_scale_factors.yaml",
                                                  f"{localdir}/params/jet_calibration.yaml",
                                                  f"{localdir}/params/skim_params.yaml",
                                                  f"{localdir}/params/sample_params.yaml",
                                                  f"{localdir}/params/event_flags.yaml",
                                                  #f"{localdir}/params/jet_scale_factors.yaml",
                                                  update = True)

# Default sequence includes: JetsCalibrator, METCalibrator, ElectronsScaleCalibrator
calibrators = default_calibrators_sequence

cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [f"{localdir}/datasets/ttHSMEFT.json",
                      f"{localdir}/datasets/ttbbSMEFT.json",
                ],
            "filter": {
                "samples":  [
                            "ttHSMEFT",
                            "ttbbSMEFT",
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
                    'ttHSMEFT' : {
                        'genMatch'   : [genMatch],
                        'non_genMatch'   : [non_genMatch]
                    },
                }
            },

        workflow = ZHbbSpanetProcessor,
        workflow_options = {"parton_jet_min_dR": 0.3,
                            "parton_jet_min_dR_postfsr": 1.0,
                            "spanet_model" : spanet_model_path,
                            "jec_pt_variation" : "",
                            "jer_variation" : "", # define the variation here
                            "dump_columns_as_arrays_per_chunk": "/cms/data/store/user/jsamudio/eftchunks" },

        skim = [precut],
        #preselections = [get_JetVetoMap(), event_selection, btag_mask, vetoE, vetoMu],
        preselections = [event_selection, btag_mask, vetoE, vetoMu],
        categories = {
            "btag_mask": [passthrough],
            },

        weights_classes = common_weights,
        weights = {
            "common": {
                "inclusive": ["genWeight"],
            }
        },
        variations = {
            "weights": {
                "common": {
                    "inclusive": [],
                }
            },
            "shape": {
                "common":{
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
            #"mAK8" : HistConf([Axis(coll="FatJetGood", field="particleNet_mass", bins = 100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
            #"zhbbtag" : HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
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
                "inclusive": [ColOut("events", outvars.common_vars+outvars.weight_vars+['sig_score', 'ttbb_score', 'ttlf_score']),
                        #ColOut(
                        #    "Parton",
                        #    ["pt", "eta", "phi", "mass", "pdgId", "provenance"]
                        #),
                        #ColOut(
                        #    "PartonMatched",
                        #    ["pt", "eta", "phi","mass", "pdgId", "provenance", "dRMatchedJet"],
                        #),
                        #ColOut(
                        #    "JetGood",
                        #    ["pt", "eta", "phi", "hadronFlavour", "btagDeepFlavB", "btag_L", "btag_M", "btag_H"],
                        #),
                        #ColOut(
                        #    "JetGoodMatched",
                        #    ["pt", "eta", "phi", "hadronFlavour", "btagDeepFlavB", "btag_L", "btag_M", "btag_H", "dRMatchedJet"],
                        #),
                        #ColOut("LeptonGood",
                        #       ["pt","eta","phi", "pdgId", "charge", "mvaTTH"],
                        #       pos_end=1, store_size=False),
                        #ColOut("MET", ["phi","pt","significance"]),
                        #ColOut("Generator",["x1","x2","id1","id2","xpdf1","xpdf2"]),
                        #ColOut("LeptonParton",["pt","eta","phi","mass","pdgId"]),
                        #ColOut("FatJetSorted",["pt", "eta", "phi", "mass", "particleNetMD_Xbb"], pos_end=1, store_size=False),
                        #ColOut("spanet_outputZ", ["ttzbb"], flatten=False),
                        #ColOut("spanet_outputH", ["ttzbb", "tthbb", "ttbb", "ttlf"], flatten=False),
                        ColOut("spanet_output", ["signal"], flatten=False), #outputs the signal from spanet
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
                #,
                #"TTZToBB": {
                #    "bycategory": {
                #        "btag_mask": [
                #            ColOut("HiggsGen",
                #                   ["pt", "eta", "phi", "mass", "pdgId"], pos_end=1, store_size=False),
                #            ]
                #        }
                #    },
                #"ttHTobb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                "ttHSMEFT": {"inclusive": [ColOut("events", outvars.sig_vars),
                            ColOut("events",["EFTfitCoefficients", "EFTfitCoefficientIndex1", "EFTfitCoefficientIndex2", "WCnames"])]
                            },
                "ttbbSMEFT": {"inclusive": [ColOut("events", outvars.bkg_vars),
                            ColOut("events",["EFTfitCoefficients", "EFTfitCoefficientIndex1", "EFTfitCoefficientIndex2", "WCnames"])]
                            },
                #"ttHToNonbb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                #"TTZToQQ": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                #"TTZToLLNuNu": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                #"TTZToBB": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                #"TTbb_Hadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTbb_SemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTbb_2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTToHadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTTo2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTToSemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"QCD_HT": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"WJetsToLNu_HT": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"DYJetsToLL_HT": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
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
