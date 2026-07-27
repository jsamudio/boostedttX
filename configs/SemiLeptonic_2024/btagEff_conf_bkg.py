from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_HLTsel_custom, get_JetVetoMap
from pocket_coffea.lib.weights.common.common import common_weights
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.parameters.histograms import *

# 1. IMPORT THE NEW PROCESSOR
from btagEffWorkflow import BTagEffProcessor 

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(custom_cut_functions)
from custom_cut_functions import *

import os
localdir = os.path.dirname(os.path.abspath(__file__))

from pocket_coffea.lib.calibrators.common.common import default_calibrators_sequence
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
                                                  update = True)

# 2. FORCE NOMINAL SHAPE CALIBRATIONS
calibrators = default_calibrators_sequence
active_shape_calibrators = ["jet_calibration"]
target_var_list = []  # Explicitly empty to prevent variations

if "jets_calibration" in parameters and "variations" in parameters.jets_calibration:
    for jet_type in parameters.jets_calibration.variations.keys():
        for year in parameters.jets_calibration.variations[jet_type].keys():
            parameters.jets_calibration.variations[jet_type][year] = target_var_list

cfg = Configurator(
        parameters = parameters,
        datasets = {
            # Keep all your JSONs here
            "jsons": [
                      f"{localdir}/datasets/ttHTobb.json",
                      f"{localdir}/datasets/TTLL.json",
                      f"{localdir}/datasets/TTNuNu.json",
                      f"{localdir}/datasets/ttHToNonbb.json",
                      f"{localdir}/datasets/TTZToQQ.json",
                      f"{localdir}/datasets/TTbb_Hadronic.json",
                      f"{localdir}/datasets/TTbb_2L2Nu.json",
                      f"{localdir}/datasets/TTbb_SemiLeptonic.json",
                      f"{localdir}/datasets/TTToHadronic.json",
                      f"{localdir}/datasets/TTTo2L2Nu.json",
                      f"{localdir}/datasets/TTToSemiLeptonic.json"
                ],
            "filter": {
                "samples":  [
                            # Process only your MC samples for efficiencies
                            # "ttHTobb",
                            # "TTLL",
                            # "TTNuNu",
                            # "ttHToNonbb",
                            # "TTZToQQ",
                            "TTbb_Hadronic",
                            "TTbb_SemiLeptonic",
                            "TTbb_2L2Nu",
                            "TTToHadronic",
                            "TTTo2L2Nu",
                            "TTToSemiLeptonic",
                ],
                "samples_exclude": [],
                "year": ['2024']
                },
                "subsamples": {
                    # 'TTbb_SemiLeptonic': {
                    #     'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    # },
                    # 'TTbb_Hadronic': {
                    #     'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    # },
                    # 'TTbb_2L2Nu': {
                    #     'tt+B'    : [get_genTtbarId_100_eq([51, 52, 53, 54, 55, 56])],
                    # },
                    # 'TTToSemiLeptonic': {
                    #     'tt+LF'   : [get_genTtbarId_100_eq(0)],
                    #     'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                    # },
                    # 'TTToHadronic': {
                    #     'tt+LF'   : [get_genTtbarId_100_eq(0)],
                    #     'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                    # },
                    # 'TTTo2L2Nu': {
                    #     'tt+LF'   : [get_genTtbarId_100_eq(0)],
                    #     'tt+C'    : [get_genTtbarId_100_eq([41, 42, 43, 44, 45, 46])],
                    # },
                #     'TTZToQQ' : {
                #         'genMatch'   : [genMatch],
                #         'non_genMatch'   : [non_genMatch]
                #     },
                #     'TTLL' : {
                #         'genMatch'   : [genMatch],
                #         'non_genMatch'   : [non_genMatch]
                #     },
                #     'ttHTobb' : {
                #         'genMatch'   : [genMatch],
                #         'non_genMatch'   : [non_genMatch]
                #     },
                #     'ttHToNonbb' : {
                #         'genMatch'   : [genMatch],
                #         'non_genMatch'   : [non_genMatch]
                #     },
                }
            },

        workflow = BTagEffProcessor,
        workflow_options = {}, 

        skim = [precut,
               get_HLTsel_custom(['HLT_Ele30_WPTight_Gsf', 'HLT_Ele115_CaloIdVT_GsfTrkIdT', 'HLT_IsoMu24', 'HLT_Mu50', 'HLT_HighPtTkMu100', 'HLT_CascadeMu100'])],
        
        # 3. CRITICAL: REMOVE `btag_mask` FROM PRESELECTIONS
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu],
        
        categories = {
            "baseline": [zh_event_cuts],
            },

        weights_classes = common_weights,
        weights = {
            "common": {
                "inclusive": ["genWeight"],
            }
        },
        variations = {
            "weights": {"common": {"inclusive": []}},
            "shape": {"common": {"inclusive": active_shape_calibrators}}
        },
        
        # 4. EMPTY VARIABLES AND COLUMNS
        variables = {
            "eff_all": HistConf(
                [
                    Axis(coll="JetGood", field="pt", bins=[20, 30, 50, 70, 100, 140, 200, 300, 600, 1000], type="variable", label="pT"),
                    Axis(coll="JetGood", field="abseta", bins=5, start=0.0, stop=2.5, type="regular", label="abseta"),
                    Axis(coll="JetGood", field="clean_flavor", bins=[-0.5, 0.5, 3.5, 4.5, 5.5], type="variable", label="flavor")
                ]
            ),
            "eff_pass": HistConf(
                [
                    Axis(coll="JetGood_Pass", field="pt", bins=[20, 30, 50, 70, 100, 140, 200, 300, 600, 1000], type="variable", label="pT"),
                    Axis(coll="JetGood_Pass", field="abseta", bins=5, start=0.0, stop=2.5, type="regular", label="abseta"),
                    Axis(coll="JetGood_Pass", field="clean_flavor", bins=[-0.5, 0.5, 3.5, 4.5, 5.5], type="variable", label="flavor")
                ]
            )
        },
        columns = {
            "common": {"inclusive": [], "bycategory": {}},
            "bysample": {}
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
        "mem_per_worker" : "4GB",
        "disk_per_worker" : "1GB",
        "exclusive"      : False,
        "chunk"          : 400000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False
    }