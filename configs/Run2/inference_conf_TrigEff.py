from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel_custom, get_JetVetoMap, eventFlags, goldenJson
from pocket_coffea.parameters.histograms import *
import spanet_inference_trigEff
from spanet_inference_trigEff import ZHbbSpanetProcessor
from pocket_coffea.lib.weights.common.common import common_weights
import outvars
from pocket_coffea.parameters.cuts import passthrough

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(spanet_inference_trigEff)
cloudpickle.register_pickle_by_value(custom_cut_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))
print(localdir)

# Import default calibrators for dynamic injection
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))  # configs/common
from met_calibration import CALIBRATORS   # default sequence with type-1 PF MET for Run 2

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
calibrators = CALIBRATORS

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


# =======================================================
# SPANET MULTI-MODEL DYNAMIC COLUMN GENERATION
# =======================================================
spanet_models_dict = {
    "baseline": {"path": "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/run2UL.onnx"},
    "overlap_NoXbb": {"path": "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/run2UL_NoXbb.onnx"},
    "overlapRemoval_NoXbb": {"path": "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/run2UL_OverlapRemovalNoXbb.onnx"},
    "overlapRemoval_FullXbb": {"path": "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/run2UL_OverlapRemovalFullXbb.onnx"},
    "overlapRemoval_MultiWP": {"path": "/cms/data/jsamudio/boosted/boostedttX/configs/SPANet/run2UL_OverlapRemovalMultiWPXbb.onnx"}
}

spanet_event_vars = []
spanet_colouts = []

for model_name in spanet_models_dict.keys():
    # 1. Grab the top-level score arrays
    spanet_event_vars.extend([
        f"sig_score_{model_name}",
        f"ttbb_score_{model_name}",
        f"ttlf_score_{model_name}"
    ])
    
    # 2. Grab the 32 encoder features
    spanet_event_vars.extend([f"encoderFeature{i+1}_{model_name}" for i in range(32)])
    
    # 3. Create the nested collection outputs (replaces your old spanet_output)
    spanet_colouts.append(ColOut(f"spanet_output_{model_name}", ["signal"], flatten=False))
# =======================================================

import awkward as ak
from pocket_coffea.lib.cut_definition import Cut

# 1. Your dictionary of triggers
hlt_paths = {
    "2016_PreVFP": [
        'HLT_Ele27_WPTight_Gsf', 'HLT_Photon175', 'HLT_Ele115_CaloIdVT_GsfTrkIdT', 
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165', 'HLT_IsoMu24', 'HLT_IsoTkMu24', 
        'HLT_Mu50', 'HLT_TkMu50'
    ],
    "2016_PostVFP": [
        'HLT_Ele27_WPTight_Gsf', 'HLT_Photon175', 'HLT_Ele115_CaloIdVT_GsfTrkIdT', 
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165', 'HLT_IsoMu24', 'HLT_IsoTkMu24', 
        'HLT_Mu50', 'HLT_TkMu50'
    ],
    "2017": [
        'HLT_Ele35_WPTight_Gsf', 'HLT_Ele32_WPTight_Gsf_L1DoubleEG', 'HLT_Photon200', 
        'HLT_Ele115_CaloIdVT_GsfTrkIdT', 'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165',
        'HLT_IsoMu27', 'HLT_Mu50', 'HLT_OldMu100', 'HLT_TkMu100'
    ],
    "2018": [
        'HLT_Ele32_WPTight_Gsf', 'HLT_Photon200', 'HLT_Ele115_CaloIdVT_GsfTrkIdT', 
        'HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165', 'HLT_IsoMu24', 'HLT_Mu50', 
        'HLT_OldMu100', 'HLT_TkMu100'
    ]
}

# 2. Extract the unique union of all trigger paths across all years (without the HLT_ prefix)
all_hlt_paths = sorted(list(set([p.replace("HLT_", "") for paths in hlt_paths.values() for p in paths])))

cfg = Configurator(
        parameters = parameters,
        datasets = {
            "jsons": [f"{localdir}/datasets/ttHTobb.json",
                      f"{localdir}/datasets/ttHToNonbb.json",
                      f"{localdir}/datasets/TTZToBB.json",
                      f"{localdir}/datasets/TTZToQQ.json",
                      f"{localdir}/datasets/TTZToLLNuNu.json",
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
                            #"ttHTobb",
                            # "ttHToNonbb",
                            # #"TTZToQQ", No scale weights?
                            # "TTZToLLNuNu",
                            # "TTZToBB",
                            #"TTbb_Hadronic",
                            "TTbb_SemiLeptonic",
                            "TTbb_2L2Nu",
                            #"TTToHadronic",
                            "TTTo2L2Nu",
                            "TTToSemiLeptonic",
                            # "WJetsToLNu",
                            # "DYJetsToLL",
                            # "SingleTop",
                            # "ttX",
                            #"DATA_EGamma",
                            #"DATA_Muon"
                            
                ],
                "samples_exclude": [],
                "year": ['2018'] # <--- FULL RUN 2 UL ENABLED
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
        calibrators = CALIBRATORS,
        workflow_options = {"parton_jet_min_dR": 0.3,
                            "parton_jet_min_dR_postfsr": 1.0,
                            "jec_pt_variation" : "",
                            "jer_variation" : "",
                            "pnet_mass_json": "/cms/data/jsamudio/boosted/boostedttX/configs/Run2/pnet_mass_jerc.json",
                            # --- MULTI-MODEL DICTIONARY INTEGRATED HERE ---
                            "spanet_models": {
                                "baseline": {
                                    "path": spanet_models_dict["baseline"]["path"],
                                    "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
                                    "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "bbtag"]
                                },
                                "overlap_NoXbb": {
                                    "path": spanet_models_dict["overlap_NoXbb"]["path"],
                                    "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
                                    "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi"]
                                },
                                "overlapRemoval_NoXbb": {
                                    "path": spanet_models_dict["overlapRemoval_NoXbb"]["path"],
                                    "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
                                    "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi"]
                                },
                                "overlapRemoval_FullXbb": {
                                    "path": spanet_models_dict["overlapRemoval_FullXbb"]["path"],
                                    "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
                                    "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "bbtag"]
                                },
                                "overlapRemoval_MultiWP": {
                                    "path": spanet_models_dict["overlapRemoval_MultiWP"]["path"],
                                    "jet_features": ["pt", "eta", "sin_phi", "cos_phi", "btagB"],
                                    "fatjet_features": ["pt", "eta", "sin_phi", "cos_phi", "xbb_wp_L", "xbb_wp_M", "xbb_wp_T"] 
                                },
                                "hlt_paths_dict": hlt_paths,
                                "all_hlt_paths": all_hlt_paths
                            }
        },
        
        skim = [eventFlags, goldenJson, precut],
        preselections = [get_JetVetoMap(), trig_event_selection_run2, vetoE, vetoMu, HEM_veto],
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
        },
        columns = {
            "common": {
                # Updated ColOut for new variables and SPANet output mapping
                "inclusive": [
                        ColOut("events", outvars.val_vars_trig+outvars.weight_vars),
                        ColOut("HLT", outvars.trig_vars2018)
                ],
                
                "bycategory": {}
            },
            "bysample": {
                # "ttHTobb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "ttHToNonbb": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToQQ": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToLLNuNu": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                # "TTZToBB": {"inclusive": [ColOut("events", outvars.sig_vars)]},
                #"TTbb_Hadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                "TTbb_SemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                "TTbb_2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                #"TTToHadronic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                "TTTo2L2Nu": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
                "TTToSemiLeptonic": {"inclusive": [ColOut("events", outvars.bkg_vars)]},
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