from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import get_nObj_min, get_HLTsel
from pocket_coffea.parameters.histograms import *
import spanet_inference_fullHad
from spanet_inference_fullHad import ZHbbSpanetProcessor
import outvars

from pocket_coffea.lib.weights.common.common import (
    common_weights,
    SF_L1prefiring,
)

import cloudpickle
import custom_cut_functions
cloudpickle.register_pickle_by_value(spanet_inference_fullHad)
cloudpickle.register_pickle_by_value(custom_cut_functions)

from custom_cut_functions import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

spanet_model_path = "/home/moanwar/ttx_hadronic/SPANet/spanet.onnx"

from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

parameters = defaults.merge_parameters_from_files(
    default_parameters,
    f"{localdir}/params/object_preselection.yaml",
    f"{localdir}/params/lepton_scale_factors.yaml",
    f"{localdir}/params/skim_params.yaml",
    f"{localdir}/params/sample_params.yaml",
    f"{localdir}/params/event_flags.yaml",
    f"{localdir}/params/triggers.yaml",
    update=True
)

# ─────────────────────────────────────────────
# Full dataset list — filter is set per run config
# ─────────────────────────────────────────────
ALL_DATASETS = [
    # Signal
    f"{localdir}/datasets/ttHTobb.json",
    f"{localdir}/datasets/ttHToNonbb.json",
    # ttZ
    f"{localdir}/datasets/TTZToBB.json",
    f"{localdir}/datasets/TTZToQQ.json",
    # ttW
    f"{localdir}/datasets/TTWJetsToLNu.json",
    f"{localdir}/datasets/TTWJetsToQQ.json",
    # TTTT, TTGJets
    f"{localdir}/datasets/TTTT.json",
    f"{localdir}/datasets/TTGJets.json",
    # tt+bb
    f"{localdir}/datasets/TTbb_Hadronic.json",
    f"{localdir}/datasets/TTbb_SemiLeptonic.json",
    f"{localdir}/datasets/TTbb_2L2Nu.json",
    # tt inclusive
    f"{localdir}/datasets/TTToHadronic.json",
    f"{localdir}/datasets/TTToSemiLeptonic.json",
    f"{localdir}/datasets/TTTo2L2Nu.json",
    # Single top
    f"{localdir}/datasets/ST_tW_top.json",
    f"{localdir}/datasets/ST_tW_antitop.json",
    f"{localdir}/datasets/ST_s_channel.json",
    f"{localdir}/datasets/ST_t_channel_top.json",
    f"{localdir}/datasets/ST_t_channel_antitop.json",
    # QCD HT
    f"{localdir}/datasets/QCD_HT500to700.json",
    f"{localdir}/datasets/QCD_HT700to1000.json",
    f"{localdir}/datasets/QCD_HT1000to1500.json",
    f"{localdir}/datasets/QCD_HT1500to2000.json",
    f"{localdir}/datasets/QCD_HT2000toInf.json",
    # W+jets
    f"{localdir}/datasets/WJetsToLNu_HT-70To100.json",
    f"{localdir}/datasets/WJetsToLNu_HT-100To200.json",
    f"{localdir}/datasets/WJetsToLNu_HT-200To400.json",
    f"{localdir}/datasets/WJetsToLNu_HT-400To600.json",
    f"{localdir}/datasets/WJetsToLNu_HT-600To800.json",
    f"{localdir}/datasets/WJetsToLNu_HT-800To1200.json",
    f"{localdir}/datasets/WJetsToLNu_HT-1200To2500.json",
    f"{localdir}/datasets/WJetsToLNu_HT-2500ToInf.json",
    f"{localdir}/datasets/WJetsToQQ_HT-200to400.json",
    f"{localdir}/datasets/WJetsToQQ_HT-400to600.json",
    f"{localdir}/datasets/WJetsToQQ_HT-600to800.json",
    f"{localdir}/datasets/WJetsToQQ_HT-800toInf.json",
    # DY+jets
    f"{localdir}/datasets/DYJetsToLL_HT-70to100.json",
    f"{localdir}/datasets/DYJetsToLL_HT-100to200.json",
    f"{localdir}/datasets/DYJetsToLL_HT-200to400.json",
    f"{localdir}/datasets/DYJetsToLL_HT-400to600.json",
    f"{localdir}/datasets/DYJetsToLL_HT-600to800.json",
    f"{localdir}/datasets/DYJetsToLL_HT-800to1200.json",
    f"{localdir}/datasets/DYJetsToLL_HT-1200to2500.json",
    f"{localdir}/datasets/DYJetsToLL_HT-2500toInf.json",
    # Z+jets
    f"{localdir}/datasets/ZJetsToQQ_HT-200to400.json",
    f"{localdir}/datasets/ZJetsToQQ_HT-400to600.json",
    f"{localdir}/datasets/ZJetsToQQ_HT-600to800.json",
    f"{localdir}/datasets/ZJetsToQQ_HT-800toInf.json",
    # Data
    f"{localdir}/datasets/JetHT.json",
    f"{localdir}/datasets/BTagCSV.json",
]

ALL_MC_SAMPLES = [
    "ttHTobb", "ttHToNonbb",
    "TTZToBB", "TTZToQQ",
    "TTWJetsToLNu", "TTWJetsToQQ",
    "TTTT", "TTGJets",
    "TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu",
    "TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu",
    "ST_tW_top", "ST_tW_antitop", "ST_s_channel",
    "ST_t_channel_top", "ST_t_channel_antitop",
    "QCD_HT",
    "WJetsToLNu_HT", "WJetsToQQ_HT",
    "DYJetsToLL_HT", "ZJetsToQQ_HT",
]

ALL_DATA_SAMPLES = ["JetHT", "BTagCSV"]

ALL_SAMPLES = ALL_MC_SAMPLES + ALL_DATA_SAMPLES

ALL_YEARS = ["2016APV", "2016", "2017", "2018"]

# ─────────────────────────────────────────────
# Subsamples
# ─────────────────────────────────────────────
SUBSAMPLES = {
    'TTbb_SemiLeptonic': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTbb_Hadronic': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTbb_2L2Nu': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTToSemiLeptonic': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTToHadronic': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTTo2L2Nu': {
        'tt+LF': [get_genTtbarId_100_eq(0)],
        'tt+C':  [get_genTtbarId_100_eq([41,42,43,44,45,46])],
        'tt+B':  [get_genTtbarId_100_eq([51,52,53,54,55,56])],
    },
    'TTZToBB':    {'genMatch': [genMatch], 'non_genMatch': [non_genMatch]},
    'TTZToQQ':    {'genMatch': [genMatch], 'non_genMatch': [non_genMatch]},
    'ttHTobb':    {'genMatch': [genMatch], 'non_genMatch': [non_genMatch]},
    'ttHToNonbb': {'genMatch': [genMatch], 'non_genMatch': [non_genMatch]},
}

# ─────────────────────────────────────────────
# Base cfg — datasets filter is intentionally empty here;
# each run config sets its own samples + years
# ─────────────────────────────────────────────
cfg = Configurator(
    parameters=parameters,
    weights_classes=common_weights + [SF_L1prefiring],
    datasets={
        "jsons": ALL_DATASETS,
        "filter": {
            "samples": ALL_SAMPLES,      # overridden per run config
            "samples_exclude": [],
            "year": ALL_YEARS,           # overridden per run config
        },
        "subsamples": SUBSAMPLES,
    },

    workflow=ZHbbSpanetProcessor,
    workflow_options={
        "parton_jet_min_dR": 0.3,
        "parton_jet_min_dR_postfsr": 1.0,
        "spanet_model": spanet_model_path,
    },

    skim=[precut, get_HLTsel(primaryDatasets=["JetHT", "BTagCSV"])],
    preselections=[event_selection],
    categories={
        "btag_mask": [btag_mask, vetoE, vetoMu],
    },

    weights={
        "common": {
            "inclusive": [
                "genWeight", "lumi", "XS", "pileup",
                "sf_L1prefiring",
                "sf_btag",
                "sf_jet_puId",
                "sf_partonshower_isr",
                "sf_partonshower_fsr",
            ]
        },
        "bysample": {
            "JetHT":   {"inclusive": []},
            "BTagCSV": {"inclusive": []},
        }
    },

    variations={
        "weights": {
            "common": {
                "inclusive": [
                    "pileup",
                    "sf_L1prefiring",
                    "sf_jet_puId",
                    "sf_btag",
                    "sf_partonshower_isr",
                    "sf_partonshower_fsr",
                ]
            }
        }
    },

    variables={
        "mAK8": HistConf([Axis(coll="FatJetGood", field="particleNet_mass",
                               bins=100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
        "zhbbtag": HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb",
                                  bins=40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
        "outZH_b1_pt": HistConf([Axis(coll="events", field="outZH_b1_pt",
                                      bins=100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)]),
    },

    columns={
        "common": {
            "inclusive": [
                ColOut("events", ["n_b_inZH", "nJetGood", "nbJetGood", "HT",
                                  "ZH_pt", "ZH_M"] + outvars.weight_vars),
                ColOut("JetGood", [
                    "pt", "eta", "phi", "mass",
                    "btagDeepFlavB", "btagDeepFlavCvL", "btagDeepFlavCvB",
                    "btagDeepFlavQG", "btag_L", "btag_M", "btag_H",
                ]),
                ColOut("MET", ["phi", "pt", "significance"]),
                ColOut("FatJetSorted", [
                    "pt", "eta", "phi", "mass",
                    "particleNetMD_Xbb", "deepTagMD_ZHbbvsQCD",
                    "deepTagMD_TvsQCD", "particleNet_TvsQCD",
                ], pos_end=1, store_size=False),
                ColOut("spanet_output", [
                    "tthbb", "ttzbb", "ttbb", "ttc", "ttlf", "qcd", "signal"
                ], flatten=False),
                ColOut("events", [
                    "spanet_t1_mass", "spanet_t2_mass",
                    "spanet_t1_w_mass", "spanet_t2_w_mass", "spanet_tt_mass"
                ]),
            ],
            "bycategory": {}
        },
        "bysample": {}
    }
)

# ─────────────────────────────────────────────
# Run options (Kodiak local)
# ─────────────────────────────────────────────
run_options = {
    "executor":        "dask/lxplus",
    "env":             "myenv",
    "cores":           4,
    "workers":         1,
    "scaleout":        50,
    "worker_image":    "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
    "queue":           "microcentury",
    "walltime":        "00:40:00",
    "mem_per_worker":  "4GB",
    "disk_per_worker": "1GB",
    "exclusive":       False,
    "chunk":           800000,
    "retries":         50,
    "treereduction":   20,
    "adapt":           False,
}
