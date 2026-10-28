"""
Shared pocket-coffea configuration of the Run 3 (2024) analysis.

Every inference_conf_<group>.py is a two-line wrapper around build_config(<group>); the groups
(which samples, data/MC, selection) are defined in samples.GROUPS and the SPANet models in
spanet_models.py.

    pocket-coffea run --cfg inference_conf_signal.py -o outputs/nominal/signal

Environment overrides:
    PC_ERAS=2024                 eras to process (default: DEFAULT_ERAS below)
    PC_SAMPLES=ttHTobb,TTLL      only these samples of the group (quick tests)
    PC_SHAPE_CALIBRATOR=...      shape calibrator to run (default: jet_calibration)
    PC_JEC_VARIATION=JES_Total   jet variation to run on top of nominal (default: nominal only)
    PC_DCTR=1                    DCTR hdamp weights dctr_hdamp_up/down for DCTR_HDAMP['samples'] (dctr.py)

Run 3 object definitions: PUPPI MET (pocket-coffea's default MET calibrator: type-1 PUPPI MET
recomputed from RawPuppiMET with the recommended JECs, data and MC), UParT AK4 b tagging, GloParT
XbbVsQCD for the fatjets and the GloParT 2-prong regressed fatjet mass.
"""
import os
import sys

import cloudpickle
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.columns_manager import ColOut
from pocket_coffea.lib.cut_functions import eventFlags, goldenJson, get_JetVetoMap
from pocket_coffea.lib.hist_manager import Axis, HistConf
from pocket_coffea.lib.weights.common.common import common_weights
from pocket_coffea.parameters import defaults
from pocket_coffea.parameters.cuts import passthrough

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))  # configs/common
import analysis_years
import custom_cut_functions
import samples
import event_variables
import spanet_inference
import spanet_models
import spanet_training
import spanet_training_workflow
import parton_matching
import xbb_wp_workflow
import btag_eff_workflow
import outvars
import dctr
from analysis_years import ERAS, BTAG_SF_WPS
from custom_cut_functions import (precut, hlt_cut, hlt_paths, event_selection, trig_event_selection, vetoE, vetoMu,
                                  btag_mask, genMatch, non_genMatch, get_genTtbarId_100_eq, get_spanet_selection, get_flag_cut)
from spanet_inference import ZHbbRun3Processor
from spanet_training_workflow import ZHbbRun3TrainingProcessor
from xbb_wp_workflow import ZHbbRun3XbbWPProcessor
from btag_eff_workflow import ZHbbRun3BtagEffProcessor

for _module in (spanet_inference, custom_cut_functions, analysis_years, samples, spanet_models,
                spanet_training, spanet_training_workflow, parton_matching, event_variables, xbb_wp_workflow, dctr,
                btag_eff_workflow):
    cloudpickle.register_pickle_by_value(_module)

LOCALDIR = os.path.dirname(os.path.abspath(__file__))

# Eras processed when PC_ERAS is not set
DEFAULT_ERAS = ['2024']

# DCTR hdamp weights (PC_DCTR=1): the ttbb has no dedicated hdamp samples in 2024 (the 5FS ttbar has)
DCTR_DIR = "/cms/data/jsamudio/boosted/boostedttX/configs/SemiLeptonic_2024/dctr_models"
DCTR_HDAMP = {
    "samples": ["TTbb_Hadronic", "TTbb_SemiLeptonic", "TTbb_2L2Nu"],
    "models": {"up": f"{DCTR_DIR}/mymodel12_hdamp_up_13.6TeV.onnx",
               "down": f"{DCTR_DIR}/mymodel12_hdamp_down_13.6TeV.onnx"},
}

# params/<name>.yaml merged on top of the pocket-coffea defaults (overrides only)
PARAMETER_FILES = ["object_preselection", "lepton_scale_factors", "jets_calibration", "lumi"]

# Jet types whose JES/JER variations are run (AK4CorrT1METJetPuppi follows AK4PFPuppi for the MET)
VARIED_JET_TYPES = ["AK4PFPuppi", "AK8PFPuppi"]

# Trigger SFs per era: derived files (TrigEff.ipynb: EleTrigSF / MuTrigSF, path relative to this
# directory or absolute), "pog" (MUO reference SF, muons only) or None (= 1).
# Derivation (orthogonal method, as in Run 2):
#   1. run the trigeff group with these settings, electron cell of TrigEff.ipynb -> ele_trig_sf_2024.json.gz
#   2. set "ele" to that file, rerun trigeff, muon cell -> muon_trig_sf_2024.json.gz
#   3. set "mu" to that file for the analysis
TRIGGER_SFS = {
    "2024": {"ele": None, "mu": "pog"},
}

# b-tag SFs per era: fixed-WP SFs (BTV UParTAK4_comb / _light, multi-WP method 1a over
# analysis_years.BTAG_SF_WPS) with our MC efficiency maps (path relative to this directory or
# absolute); None = btag_sf 1. Derivation:
#   1. run_pipeline.py inference --groups btageff
#   2. python3 derive_btag_eff.py -> btag_eff_2024.json.gz; set "eff" to that file
BTAG_SFS = {
    "2024": {"eff": None},
}

# Efficiency-map histogram of the btageff group (AK4 jets of the preselection, btag_eff_workflow.py)
BTAG_EFF_PT_BINS = [30, 50, 70, 100, 140, 200, 300, 600, 1000]
BTAG_EFF_ABSETA_BINS = [0.0, 0.6, 1.2, 1.8, 2.4]
BTAG_EFF_HISTOGRAMS = {
    "btag_eff": HistConf([
        Axis(coll="JetGood", field="hadronFlavour", type="intcat", bins=[0, 4, 5], name="flavor",
             label="Hadron flavour"),
        Axis(coll="JetGood", field="pt", type="variable", bins=BTAG_EFF_PT_BINS, name="pt",
             label=r"$p_T$ [GeV]"),
        Axis(coll="JetGood", field="abseta", type="variable", bins=BTAG_EFF_ABSETA_BINS, name="abseta",
             label=r"$|\eta|$"),
        Axis(coll="JetGood", field="btag_level", type="int", start=0, stop=len(BTAG_SF_WPS) + 1,
             name="level", label=f"b-tag level (0 = untagged, 1..{len(BTAG_SF_WPS)} = {'/'.join(BTAG_SF_WPS)})"),
    ]),
}

HISTOGRAMS = {
    "mAK8": HistConf([Axis(coll="FatJetGood", field="mass", bins=100, start=0, stop=200,
                           label=r"$M_{GloParT}$ [GeV]")]),
    "zhbbtag_sorted": HistConf([Axis(coll="FatJetSorted", field="btagBB", bins=40, start=0, stop=1,
                                     label=r"GloParT $X_{bb}$ vs QCD", pos=0)]),
}


def _resolve_trigger_sfs(eras):
    out = {}
    for era in eras:
        out[era] = {}
        for flavour, value in TRIGGER_SFS.get(era, {}).items():
            if value not in (None, "pog") and not os.path.isabs(value):
                value = os.path.join(LOCALDIR, value)
            if value not in (None, "pog") and not os.path.isfile(value):
                raise FileNotFoundError(f"Trigger SF file for {flavour} {era} not found: {value}")
            out[era][flavour] = value
    return out


def _resolve_btag_sfs(eras):
    """Workflow option 'btag_sfs': the working points and the efficiency-map file per era."""
    eff = {}
    for era in eras:
        value = BTAG_SFS.get(era, {}).get("eff")
        if value and not os.path.isabs(value):
            value = os.path.join(LOCALDIR, value)
        if value and not os.path.isfile(value):
            raise FileNotFoundError(f"b-tag efficiency maps for {era} not found: {value}")
        eff[era] = value
    return {"wps": list(BTAG_SF_WPS), "eff": eff}


def _check_btag_wps(parameters):
    """The b-tag SFs must cover every working point whose information the analysis uses."""
    used = set(spanet_training.BTAG_COUNT_WPS) | {parameters.object_preselection.btag.wp}
    missing = used - set(BTAG_SF_WPS)
    if missing:
        raise ValueError(f"b-tag working points {sorted(missing)} are used but not in "
                         f"analysis_years.BTAG_SF_WPS {BTAG_SF_WPS}: their SFs would be missing")


def _env_list(name):
    value = os.environ.get(name, "").strip()
    return [v.strip() for v in value.split(",") if v.strip()] if value else None


def available_jet_variations(parameters, era):
    """Jet variations pocket-coffea knows for an era (JES_Total, JER, and the full JES sources)."""
    default = parameters.default_jets_calibration.variations
    out = set()
    for key in ("total_variation", "full_variations"):
        node = default[key].get("AK4PFPuppi", {})
        if era in node:
            out |= set(node[era])
    return out


def dctr_options():
    """Workflow option 'dctr_hdamp' (configs/common/dctr.py): DCTR_HDAMP with PC_DCTR=1, else None (off)."""
    if not os.environ.get("PC_DCTR"):
        return None
    missing = [var for var, path in DCTR_HDAMP["models"].items() if not path]
    if missing:
        raise ValueError(f"PC_DCTR=1 but no DCTR model for {missing}: set DCTR_HDAMP['models'] in inference_common.py")
    print(f"--> DCTR hdamp weights for {DCTR_HDAMP['samples']}")
    return {"samples": list(DCTR_HDAMP["samples"]), "models": dict(DCTR_HDAMP["models"])}


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
    parameters.jets_calibration.variations = {
        jet_type: {era: list(jet_variations) for era in ERAS} for jet_type in VARIED_JET_TYPES}
    print(f"--> Shape calibrator: {shape_calibrator}, jet variations: {jet_variations or 'nominal only'}")
    return parameters, [shape_calibrator]


def build_spanet_subsamples():
    subsamples = {}
    for sample, opts in spanet_training.SAMPLES.items():
        subs = opts.get("subsamples")
        if subs == "ttbar":
            subsamples[sample] = {name: [get_genTtbarId_100_eq(ids)] for name, ids in samples.TTBAR_SPLIT.items()}
        elif subs:
            subsamples[sample] = {name: [get_flag_cut(flag)] for name, flag in subs.items()}
    return subsamples


def build_spanet_columns(category):
    """Columns of the SPANet training inputs (spanet_training.EXPORT), MC only."""
    cols = []
    for collection, spec in spanet_training.EXPORT.items():
        fields = spanet_training.exported_fields(collection)
        if spec.get("max") == 1:
            cols.append(ColOut(collection, fields, pos_end=1, store_size=False))
        else:
            cols.append(ColOut(collection, fields, pos_end=spec.get("max")))
    return {"common": {"inclusive": [], "bycategory": {category: cols}}, "bysample": {}}


# Xbb working points: all fatjets of the preselection (derive_xbb_wps.py)
XBBWP_FATJET_COLUMNS = ["pt", "eta", "phi", "mass", "btagBB", "nBHadrons", "nCHadrons", "matched", "n_boson_b"]


def build_xbbwp_subsamples():
    """tt+LF only for the ttbar samples, Z -> bb only for TTZToQQ."""
    subsamples = {s: {"tt+LF": [get_genTtbarId_100_eq(samples.TTBAR_SPLIT["tt+LF"])]}
                  for s in samples.XBBWP_MISTAG_SAMPLES}
    subsamples["TTZToQQ"] = {"Zbb": [get_flag_cut("Zbb")]}
    return subsamples


def build_xbbwp_columns(category):
    cols = [ColOut("FatJetGood", XBBWP_FATJET_COLUMNS),
            ColOut("events", ["nJetGood", "nbJetGood", "nLeptonGood"])]
    return {"common": {"inclusive": [], "bycategory": {category: cols}}, "bysample": {}}


def dataset_jsons(sample_list):
    """Dataset files: the nominal ones, and those of the variation samples that are processed
    (built from datasets/ttbar_variations_definitions.json, so they may not exist yet)."""
    names = list(samples.DATASET_JSONS)
    for name in samples.TTBAR_VARIATION_SAMPLES:
        if name in sample_list:
            if not os.path.exists(f"{LOCALDIR}/datasets/{name}.json"):
                raise FileNotFoundError(f"datasets/{name}.json not built: pocket-coffea build-datasets "
                                        f"--cfg datasets/ttbar_variations_definitions.json")
            names.append(name)
    return [f"{LOCALDIR}/datasets/{name}.json" for name in names]


def build_subsamples():
    subsamples = {}
    for sample in samples.TTBAR_SAMPLES + samples.TTBAR_VARIATION_SAMPLES:
        subsamples[sample] = {name: [get_genTtbarId_100_eq(ids)] for name, ids in samples.TTBAR_SPLIT.items()}
    for sample in samples.SIGNAL_SAMPLES:
        subsamples[sample] = {"genMatch": [genMatch], "non_genMatch": [non_genMatch]}
    return subsamples


def build_columns(sample_list, category, selection):
    """Columns common to data and MC per category; MC-only columns (weights, gen flags) per sample,
    so groups mixing data and MC (trigeff) work."""
    mc_samples = [s for s in sample_list if s not in samples.DATA_SAMPLES]
    bysample = {}

    if selection == "trigeff":
        common = [ColOut("events", outvars.trigeff_vars),
                  ColOut("HLT", sorted({p.replace("HLT_", "") for paths in hlt_paths.values() for p in paths}))]
        mc_cols = [ColOut("events", outvars.weight_vars)]
    else:
        models = spanet_models.ACTIVE_MODELS
        spanet_vars = [c for m in models for c in spanet_models.score_columns(m) + spanet_models.encoder_columns(m)]
        data_safe = [v for v in outvars.common_vars if v not in outvars.mc_only_vars]
        common = ([ColOut("events", data_safe + outvars.validation_vars + spanet_vars)]
                  + [ColOut(f"spanet_output_{m}", ["signal"], flatten=False) for m in models])
        mc_cols = [ColOut("events", outvars.mc_only_vars + outvars.weight_vars)]

    for sample in mc_samples:
        cols = list(mc_cols)
        if selection != "trigeff":
            extra = {"sig": outvars.sig_vars, "bkg": outvars.bkg_vars}.get(samples.extra_column_group(sample), [])
            if extra:
                cols.append(ColOut("events", extra))
            if sample in samples.EFT_SAMPLES:
                cols.append(ColOut("events", samples.EFT_COLUMNS))
        bysample[sample] = {"bycategory": {category: cols}}

    return {"common": {"inclusive": [], "bycategory": {category: common}}, "bysample": bysample}


def build_config(group, eras=None, sample_list=None):
    """Configurator for one sample group of samples.GROUPS."""
    opts = samples.GROUPS[group]
    eras = eras or _env_list("PC_ERAS") or DEFAULT_ERAS
    sample_list = sample_list or _env_list("PC_SAMPLES") or opts["samples"]
    unknown = [e for e in eras if e not in ERAS]
    if unknown:
        raise ValueError(f"Unknown era(s) {unknown}; choose from {ERAS}")
    selection = opts.get("selection", "analysis")
    print(f"--> Group '{group}' ({selection} selection): samples {sample_list}, eras {eras}, "
          f"SPANet models {spanet_models.ACTIVE_MODELS}")

    parameters, shape_calibrators = load_parameters(eras)
    _check_btag_wps(parameters)

    workflow = ZHbbRun3Processor
    weights = ["genWeight"]
    subsamples = build_subsamples()
    if selection == "trigeff":
        # ttbar e-mu without trigger requirement; the notebook applies reference / target triggers
        skim = [eventFlags, goldenJson, precut]
        preselections = [get_JetVetoMap(), trig_event_selection, vetoE, vetoMu]
        category = "trigeff"
        categories = {category: [passthrough]}
    elif selection == "spanet":
        # analysis preselection + the training category of spanet_training.SELECTION
        spanet_training.validate()
        data = [s for s in sample_list if s in samples.DATA_SAMPLES]
        if data:
            raise ValueError(f"The SPANet training inputs need gen information, not data: {data}")
        skim = [eventFlags, goldenJson, precut, hlt_cut]
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu]
        category = spanet_training.CATEGORY
        spanet_cut = get_spanet_selection(spanet_training.SELECTION)
        categories = {category: [spanet_cut]}
        workflow = ZHbbRun3TrainingProcessor
        weights = list(spanet_training.WEIGHTS)
        subsamples = build_spanet_subsamples()
    elif selection == "xbbwp":
        # analysis preselection, no b-tag category and no Xbb cut (nbJetGood is exported)
        data = [s for s in sample_list if s in samples.DATA_SAMPLES]
        if data:
            raise ValueError(f"The Xbb working points need gen information, not data: {data}")
        skim = [eventFlags, goldenJson, precut, hlt_cut]
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu]
        category = "presel"
        categories = {category: [passthrough]}
        workflow = ZHbbRun3XbbWPProcessor
        weights = ["genWeight", "lumi", "XS", "pileup"]
        subsamples = build_xbbwp_subsamples()
    elif selection == "btageff":
        # analysis preselection, before the b-tag requirement; histograms only (BTAG_EFF_HISTOGRAMS)
        data = [s for s in sample_list if s in samples.DATA_SAMPLES]
        if data:
            raise ValueError(f"The b-tag efficiency maps need the hadron flavour, not data: {data}")
        skim = [eventFlags, goldenJson, precut, hlt_cut]
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu]
        category = "presel"
        categories = {category: [passthrough]}
        workflow = ZHbbRun3BtagEffProcessor
        weights = ["genWeight", "lumi", "XS", "pileup"]
        subsamples = {s: {name: [get_genTtbarId_100_eq(ids)] for name, ids in samples.TTBAR_SPLIT.items()}
                      for s in samples.TTBAR_SAMPLES}
    else:
        skim = [eventFlags, goldenJson, precut, hlt_cut]
        preselections = [get_JetVetoMap(), event_selection, vetoE, vetoMu]
        category = "btag_mask"
        categories = {category: [btag_mask]}

    workflow_options = {
        "process_labels": dict(samples.PROCESS_LABELS),
        "signal_samples": samples.SIGNAL_SAMPLES + ["ttHSMEFT"],
        "trigger_sfs": _resolve_trigger_sfs(eras),
        "spanet_models": spanet_models.active_models() if selection == "analysis" else {},
        "prefilter_btag_mask": selection == "analysis",
        "dctr_hdamp": dctr_options() if selection == "analysis" else None,
        "btag_sfs": _resolve_btag_sfs(eras),
    }
    if selection == "spanet":
        workflow_options.update({
            "clean_jets": spanet_training.CLEAN_JETS,
            "matching": dict(spanet_training.MATCHING),
            "extra_b_samples": [s for s, o in spanet_training.SAMPLES.items() if o.get("extra_b_as_boson")],
            # plain dict (the workflow options are saved as JSON); the workflow builds the Cut
            "spanet_selection": dict(spanet_training.SELECTION),
            "event_features": {"btag_wps": list(spanet_training.BTAG_COUNT_WPS),
                               "count_clean": spanet_training.COUNT_CLEAN_JETS},
        })
    if opts.get("dump_dir"):
        workflow_options["dump_columns_as_arrays_per_chunk"] = opts["dump_dir"].format(era="_".join(eras))

    return Configurator(
        parameters=parameters,
        datasets={
            "jsons": dataset_jsons(sample_list),
            "filter": {"samples": sample_list, "samples_exclude": [], "year": eras},
            "subsamples": subsamples,
        },
        workflow=workflow,
        workflow_options=workflow_options,
        # MET filters and golden JSON (data) are standard pocket-coffea skim cuts
        skim=skim,
        preselections=preselections,
        categories=categories,
        weights_classes=common_weights,
        weights={"common": {"inclusive": weights}},
        variations={
            "weights": {"common": {"inclusive": []}},
            "shape": {"common": {"inclusive": shape_calibrators}},
        },
        variables=BTAG_EFF_HISTOGRAMS if selection == "btageff" else HISTOGRAMS,
        columns=(build_spanet_columns(category) if selection == "spanet"
                 else build_xbbwp_columns(category) if selection == "xbbwp"
                 else {"common": {"inclusive": [], "bycategory": {}}, "bysample": {}} if selection == "btageff"
                 else build_columns(sample_list, category, selection)),
    )
