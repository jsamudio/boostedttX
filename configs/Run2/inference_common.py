"""
Shared pocket-coffea configuration of the SPANet inference step.

Every inference_conf_<group>.py is a two-line wrapper around build_config(<group>);
the groups (which samples, data/MC, ...) are defined in samples.GROUPS and
the SPANet models in spanet_models.py.

    pocket-coffea run --cfg inference_conf_signal.py -o output_signal

Environment overrides (no need to edit any file for these):
    PC_ERAS=2017,2018            eras to process (default: DEFAULT_ERAS below)
    PC_SAMPLES=ttHTobb,TTZToBB   only these samples of the group (quick tests)
    PC_SHAPE_CALIBRATOR=...      shape calibrator to run (default: jet_calibration)
    PC_JEC_VARIATION=JES_Total   jet variation to run on top of nominal (default: nominal only)
"""
import os

import cloudpickle
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import eventFlags, goldenJson
from pocket_coffea.lib.hist_manager import Axis, HistConf
from pocket_coffea.lib.weights.common.common import common_weights
from pocket_coffea.parameters import defaults

import analysis_years
import custom_cut_functions
import samples
import spanet_inference
import spanet_models
import outvars
from analysis_years import xbb_wp_cuts, ERAS
from custom_cut_functions import (precut, hlt_cut, get_JetVetoMap, event_selection, vetoE, vetoMu,
                                  HEM_veto, btag_mask, genMatch, non_genMatch, get_genTtbarId_100_eq,
                                  trig_event_selection_run2, HLT_PATHS)
from pocket_coffea.parameters.cuts import passthrough
from spanet_inference import ZHbbSpanetProcessor
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))  # configs/common
from met_calibration import CALIBRATORS   # type-1 PF MET for Run 2 (configs/common/met_calibration.py)

for _module in (spanet_inference, custom_cut_functions, analysis_years, samples, spanet_models):
    cloudpickle.register_pickle_by_value(_module)

LOCALDIR = os.path.dirname(os.path.abspath(__file__))

# Eras processed when PC_ERAS is not set
DEFAULT_ERAS = ['2018']

# params/<name>.yaml merged on top of the pocket-coffea defaults (overrides only).
# skim_params holds the lepton selection cuts read by object_cleaning_functions.lep_sel / soft_lep_sel.
PARAMETER_FILES = ["object_preselection", "skim_params", "lepton_scale_factors", "jets_calibration",
                   "event_flags", "jet_scale_factors"]

# Jet types whose JES/JER variations are run (AK4CorrT1METJetchs follows AK4PFchs for the type-1 MET)
VARIED_JET_TYPES = ["AK4PFchs", "AK8PFPuppi"]

# Kept for reference / parity with the old configs: the pocket-coffea runner takes its run options
# from the command line or a --custom-run-options YAML, not from this dictionary.
RUN_OPTIONS = {
    "executor"       : "dask/lxplus",
    "env"            : "myenv",
    "cores"          : 4,
    "workers"        : 1,
    "scaleout"       : 50,
    "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
    "queue"          : "microcentury",
    "walltime"       : "00:40:00",
    "mem_per_worker" : "4GB",
    "disk_per_worker": "1GB",
    "exclusive"      : False,
    "chunk"          : 400000,
    "retries"        : 50,
    "treereduction"  : 20,
    "adapt"          : False,
}

HISTOGRAMS = {
    "mAK8": HistConf([Axis(coll="FatJetGood", field="particleNet_mass", bins=100, start=0, stop=200, label=r"$M_{pNet}$ [GeV]")]),
    "zhbbtag": HistConf([Axis(coll="FatJetGood", field="particleNetMD_Xbb", bins=40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
    "zhbbtag_sorted": HistConf([Axis(coll="FatJetSorted", field="particleNetMD_Xbb", bins=40, start=0, stop=1, label=r"$Xbb_{pNet}$", pos=0)]),
    "outZH_b1_pt": HistConf([Axis(coll="events", field="outZH_b1_pt", bins=100, start=0, stop=200, label=r"$Xbb_{pNet}$", pos=0)]),
}


def _env_list(name):
    value = os.environ.get(name, "").strip()
    return [v.strip() for v in value.split(",") if v.strip()] if value else None


def available_jet_variations(parameters, era):
    """Jet variations pocket-coffea knows for an era (JES_Total, JER, JES_Regrouped_*)."""
    default = parameters.default_jets_calibration.variations
    return set(default.total_variation.AK4PFchs[era]) | set(default.full_variations.AK4PFchs[era])


def load_parameters(eras):
    """Default pocket-coffea parameters + params/*.yaml, with the shape variation from the environment."""
    defaults.register_configuration_dir("config_dir", LOCALDIR + "/params")
    parameters = defaults.merge_parameters_from_files(
        defaults.get_default_parameters(),
        *[f"{LOCALDIR}/params/{name}.yaml" for name in PARAMETER_FILES],
        update=True)

    shape_calibrator = os.environ.get("PC_SHAPE_CALIBRATOR") or "jet_calibration"
    jec_var = os.environ.get("PC_JEC_VARIATION")
    jet_variations = [jec_var] if (shape_calibrator == "jet_calibration" and jec_var
                                   and jec_var.lower() != "nominal") else []
    for era in eras:
        unknown = set(jet_variations) - available_jet_variations(parameters, era)
        if unknown:
            raise ValueError(f"PC_JEC_VARIATION={jec_var} is not a jet variation for {era}; available: "
                             f"{sorted(available_jet_variations(parameters, era))}")
    # Only the requested variation (or none) is run in this job
    parameters.jets_calibration.variations = {
        jet_type: {era: list(jet_variations) for era in ERAS} for jet_type in VARIED_JET_TYPES}
    print(f"--> Shape calibrator: {shape_calibrator}, jet variations: {jet_variations or 'nominal only'}")
    return parameters, [shape_calibrator]


def build_subsamples():
    subsamples = {}
    for sample in samples.TTBAR_SAMPLES:
        subsamples[sample] = {name: [get_genTtbarId_100_eq(ids)] for name, ids in samples.TTBAR_SPLIT.items()}
    for sample in samples.SIGNAL_SAMPLES:
        subsamples[sample] = {"genMatch": [genMatch], "non_genMatch": [non_genMatch]}
    return subsamples


def build_columns(sample_list, is_data):
    models = spanet_models.ACTIVE_MODELS
    spanet_vars = [c for m in models for c in spanet_models.score_columns(m) + spanet_models.encoder_columns(m)]
    spanet_colouts = [ColOut(f"spanet_output_{m}", ["signal"], flatten=False) for m in models]

    if is_data:
        event_vars = [v for v in outvars.common_vars if v not in outvars.mc_only_vars] + outvars.validation_vars
    else:
        event_vars = outvars.common_vars + outvars.validation_vars + outvars.weight_vars

    bysample = {}
    for sample in sample_list:
        extra = {"sig": outvars.sig_vars + outvars.pdf_vars,
                 "bkg": outvars.bkg_vars + outvars.pdf_vars}.get(samples.extra_column_group(sample), [])
        cols = [ColOut("events", extra)] if extra else []
        if sample in samples.EFT_SAMPLES:
            cols.append(ColOut("events", samples.EFT_COLUMNS))
        if cols:
            bysample[sample] = {"inclusive": cols}

    return {
        "common": {
            "inclusive": [],
            "bycategory": {"btag_mask": [ColOut("events", event_vars + spanet_vars)] + spanet_colouts},
        },
        "bysample": bysample,
    }


# HLT paths of all Run 2 eras (without "HLT_"), exported by the trigger-efficiency groups;
# paths that do not exist in an era are filled with False by the processor
TRIGEFF_HLT_COLUMNS = sorted({p.replace("HLT_", "") for era in HLT_PATHS.values()
                              for paths in era.values() for p in paths})


def build_trigeff_columns(sample_list):
    """Columns read by TrigEff.ipynb: lepton kinematics, jet counts and HLT bits for all samples,
    plus the MC weights for the MC samples. The category keeps the name 'btag_mask' for the notebook."""
    common = [ColOut("events", outvars.val_vars_trig + ["n_ak4jets", "n_b_outZH"]),
              ColOut("HLT", TRIGEFF_HLT_COLUMNS)]
    mc_only = [ColOut("events", outvars.weight_vars + outvars.bkg_vars)]
    bysample = {s: {"inclusive": mc_only} for s in sample_list if s not in samples.DATA_SAMPLES}
    return {"common": {"inclusive": [], "bycategory": {"btag_mask": common}}, "bysample": bysample}


def build_config(group, eras=None, sample_list=None):
    """Configurator for one sample group of samples.GROUPS."""
    opts = samples.GROUPS[group]
    eras = eras or _env_list("PC_ERAS") or DEFAULT_ERAS
    sample_list = sample_list or _env_list("PC_SAMPLES") or opts["samples"]
    unknown = [e for e in eras if e not in ERAS]
    if unknown:
        raise ValueError(f"Unknown era(s) {unknown}; choose from {ERAS}")
    is_data = opts.get("is_data", False)
    selection = opts.get("selection", "analysis")
    print(f"--> Inference group '{group}' ({selection}): samples {sample_list}, eras {eras}, "
          f"SPANet models {spanet_models.ACTIVE_MODELS if selection == 'analysis' else 'none'}")

    parameters, shape_calibrators = load_parameters(eras)

    # HEM_veto only acts on 2018 (MC: the HEM-affected lumi fraction, data: run >= 319077)
    if selection == "trigeff":
        # ttbar e-mu, no trigger requirement (TrigEff.ipynb applies the reference / target triggers)
        skim = [eventFlags, goldenJson, precut]
        preselections = [get_JetVetoMap(), trig_event_selection_run2, vetoE, vetoMu, HEM_veto]
        categories = {"btag_mask": [passthrough]}
        columns = build_trigeff_columns(sample_list)
    else:
        skim = [eventFlags, goldenJson, precut, hlt_cut]
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu, HEM_veto]
        # The b-tag requirement is a category (not a preselection) so that the processor can
        # accumulate the b-tag SF normalization sums on the events before it.
        categories = {"btag_mask": [btag_mask]}
        columns = build_columns(sample_list, is_data)

    workflow_options = {
        "parton_jet_min_dR": 0.3,
        "parton_jet_min_dR_postfsr": 1.0,
        "pnet_mass_json": f"{LOCALDIR}/pnet_mass_jerc.json",
        "sf_dir": LOCALDIR,
        "spanet_models": spanet_models.active_models() if selection == "analysis" else {},
        "xbb_wp_cuts": {era: xbb_wp_cuts(era) for era in ERAS},
        "btag_calib": not is_data and selection == "analysis",
        "prefilter_btag_mask": selection == "analysis",
        "trigger_reference_muon": opts.get("reference_muon", False),
        "hlt_columns": TRIGEFF_HLT_COLUMNS if selection == "trigeff" else [],
    }
    if opts.get("dump_dir"):
        workflow_options["dump_columns_as_arrays_per_chunk"] = opts["dump_dir"].format(era="_".join(eras))

    return Configurator(
        parameters=parameters,
        datasets={
            "jsons": [f"{LOCALDIR}/datasets/{name}.json" for name in samples.DATASET_JSONS],
            "filter": {
                "samples": sample_list,
                "samples_exclude": [],
                "year": eras,
            },
            "subsamples": build_subsamples(),
        },
        workflow=ZHbbSpanetProcessor,
        calibrators=CALIBRATORS,
        workflow_options=workflow_options,
        # MET filters (params/event_flags.yaml) and golden JSON (data) are standard pocket-coffea skim cuts
        skim=skim,
        preselections=preselections,
        categories=categories,
        weights_classes=common_weights,
        weights={"common": {"inclusive": ["genWeight"]}},
        variations={
            "weights": {"common": {"inclusive": []}},
            "shape": {"common": {"inclusive": shape_calibrators}},
        },
        variables=HISTOGRAMS,
        columns=columns,
    )
