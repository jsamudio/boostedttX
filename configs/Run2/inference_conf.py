from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel_custom, get_JetVetoMap
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

# Import default calibrators for dynamic injection
from pocket_coffea.lib.calibrators.common.common import default_calibrators_sequence

spanet_model_path = "/cms/data/jsamudio/boosted/boostedttX/configs/spanetBalanced24.onnx"

# Load defaults
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/lepton_scale_factors.yaml",
                                                  f"{localdir}/params/jets_calibration.yaml",
                                                  f"{localdir}/params/skim_params.yaml",
                                                  f"{localdir}/params/sample_params.yaml",
                                                  f"{localdir}/params/event_flags.yaml",
                                                  f"{localdir}/params/jet_scale_factors.yaml",
                                                  update = True)

# --- DYNAMIC CALIBRATOR INJECTION ---
# Default sequence includes: JetsCalibrator, METCalibrator, ElectronsScaleCalibrator, Muon Rochester
calibrators = default_calibrators_sequence

# Grab environment variables (defaults to None if not set)
target_shape_calibrator = os.environ.get("PC_SHAPE_CALIBRATOR", None)
target_jec_var = os.environ.get("PC_JEC_VARIATION", None)

# 1. Determine which calibrator block to run
if target_shape_calibrator == "jet_calibration":
    active_shape_calibrators = ["jet_calibration"]
elif target_shape_calibrator:
    active_shape_calibrators = [target_shape_calibrator]
else:
    active_shape_calibrators = ["jet_calibration"]
    print("--> No target shape calibrator set. Defaulting to jet_calibration.")

print(f"--> Active Shape Calibrators: {active_shape_calibrators}")

# 2. Determine the Jet variations list
if target_shape_calibrator == "jet_calibration" and target_jec_var and target_jec_var.lower() != "nominal":
    target_var_list = [target_jec_var]
else:
    target_var_list = [] 

# 3. Safely update the nested dictionary structure for jets
if "jets_calibration" in parameters and "variations" in parameters.jets_calibration:
    for jet_type in parameters.jets_calibration.variations.keys():
        for year in parameters.jets_calibration.variations[jet_type].keys():
            # Update the inner list without breaking the dict structure
            parameters.jets_calibration.variations[jet_type][year] = target_var_list

print(f"--> Jet Systematic Variations set to: {target_var_list}")


cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [f"{localdir}/datasets/ttHTobb.json",
                      f"{localdir}/datasets/ttHToNonbb.json",
                      f"{localdir}/datasets/TTZToBB.json",
                      f"{localdir}/datasets/TTZToQQ.json",
                      #f"{localdir}/datasets/TTZToLLNuNu.json",
                      f"{localdir}/datasets/TTbb_Hadronic.json",
                      f"{localdir}/datasets/TTbb_2L2Nu.json",
                      f"{localdir}/datasets/TTbb_SemiLeptonic.json",
                      f"{localdir}/datasets/TTToHadronic.json",
                      f"{localdir}/datasets/TTTo2L2Nu.json",
                      f"{localdir}/datasets/TTToSemiLeptonic.json",
                      f"{localdir}/datasets/QCD.json",
                      f"{localdir}/datasets/WJetsToLNu.json",
                      f"{localdir}/datasets/DYJetsToLL.json",
                      f"{localdir}/datasets/SingleTop.json",
                      f"{localdir}/datasets/ttX.json",
                      f"{localdir}/datasets/DATA_EGamma.json",
                      f"{localdir}/datasets/DATA_Muon.json",
                      
                ],
            "filter": {
                "samples":  [
                            "ttHTobb",
                            # "ttHToNonbb",
                            # #"TTZToQQ", No scale weights?
                            # "TTZToLLNuNu",
                            # "TTZToBB",
                            # "TTbb_Hadronic",
                            # "TTbb_SemiLeptonic",
                            # "TTbb_2L2Nu",
                            # "TTToHadronic",
                            # "TTTo2L2Nu",
                            # "TTToSemiLeptonic",
                            # "WJetsToLNu",
                            # "DYJetsToLL",
                            # "SingleTop",
                            # "ttX",
                            "DATA_EGamma",
                            "DATA_Muon"
                            
                ],
                "samples_exclude": [],
                "year": ['2018'] # Updated to 2018
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

        workflow = ZHbbSpanetProcessor,
        workflow_options = {"parton_jet_min_dR": 0.3,
                            "parton_jet_min_dR_postfsr": 1.0,
                            "spanet_model" : spanet_model_path,
                            "jec_pt_variation" : "",
                            "jer_variation" : "",
                            "pnet_mass_json": "/cms/data/jsamudio/boosted/boostedttX/configs/Run2/pnet_mass_jerc.json",},

        # Updated Skim block to include the 2018 specific single-e and single-mu triggers
        skim = [
            precut,
            get_HLTsel_custom([
                'HLT_Ele32_WPTight_Gsf', 
                'HLT_Photon200', 
                'HLT_Ele115_CaloIdVT_GsfTrkIdT', 
                'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
                'HLT_IsoMu24', 
                'HLT_Mu50', 
                'HLT_OldMu100', 
                'HLT_TkMu100'
            ])
        ],
        preselections = [get_JetVetoMap(), event_selection, btag_mask, vetoE, vetoMu],
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
                    "inclusive": active_shape_calibrators, # Injecting dynamic calibrator string here
            }
        }
        },
        variables = {
            "mAK8" : HistConf([Axis(coll="FatJetGood", field="particleNet_mass", bins = 100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
            "zhbbtag" : HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            "zhbbtag_sorted" : HistConf([Axis(coll="FatJetSorted", field="particleNetMD_Xbb", bins = 40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
            "outZH_b1_pt" : HistConf([Axis(coll="events", field="outZH_b1_pt", bins = 100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)])
        },
        columns = {
            "common": {
                # Updated ColOut for new variables and SPANet output mapping
                "inclusive": [
                        ColOut("events", outvars.common_vars+outvars.validation_vars),#+['sig_score', 'ttbb_score', 'ttlf_score']),
                        ColOut("spanet_output", ["signal"], flatten=False), 
                ],
                
                "bycategory": {}
            },
            "bysample": {
                # "ttHTobb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "ttHToNonbb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToQQ": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToLLNuNu": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToBB": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTbb_Hadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                # "TTbb_SemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                # "TTbb_2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                # "TTToHadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                # "TTTo2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                # "TTToSemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
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