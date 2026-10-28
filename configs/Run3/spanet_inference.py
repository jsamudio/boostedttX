import awkward as ak
import numpy as np
import vector

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.scale_factors import sf_ele_reco, sf_ele_id, sf_pileup_reweight
from object_cleaning_functions import (soft_lep_sel, lep_sel, bjet_sel, qjet_sel, lep_softlep_combo,
                                       jet_selection_custom)
from custom_cut_functions import sortbyscore, btag_mask
from cand_helper import zh_helper
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig, match_gen_bb_to_zh
from weight_handler import calc_weight
from custom_scale_factors import (set_fatjet_mass, sf_mu, sf_partonshower_isr, sf_partonshower_fsr,
                                  sf_scaleweights, sf_ele_trig, sf_mu_trig)
from get_controlvars import get_controlvars
from spanet_models import SPANET_CLASSES
from event_variables import event_variables
from dctr import apply_dctr_hdamp, set_default_dctr_columns

vector.register_awkward()

_LOCAL_CACHE = {}


def _worker_cache():
    """Dict that survives between chunks: on the dask worker if there is one, else module level."""
    try:
        from dask.distributed import get_worker
        worker = get_worker()
    except (ImportError, ValueError):
        return _LOCAL_CACHE
    if not hasattr(worker, "zhbb_cache"):
        worker.zhbb_cache = {}
    return worker.zhbb_cache


def _onnx_session(path):
    cache = _worker_cache()
    if ("onnx", path) not in cache:
        import onnxruntime as ort
        opts = ort.SessionOptions()
        opts.intra_op_num_threads = 1
        opts.inter_op_num_threads = 1
        opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        cache[("onnx", path)] = ort.InferenceSession(path, sess_options=opts, providers=['CPUExecutionProvider'])
    return cache[("onnx", path)]


def _ones(events):
    return np.ones(len(events))


class ZHbbRun3Processor(BaseProcessorABC):
    '''
    Run 3 (2024) workflow: PUPPI MET, UParT AK4 b tagging, GloParT XbbVsQCD ranking of the fatjets
    and the GloParT 2-prong regressed fatjet mass. Per chunk and shape variation, after the
    preselection:

      stage 1 (all preselected events)
          object/event scale factors, process labels, gen matching
      prefilter to the btag_mask category (workflow option 'prefilter_btag_mask')
      stage 2
          ZH candidate variables, control variables, SPANet inference (if models are given)

    Workflow options:
      process_labels      - {sample: label} for samples not labelled from genTtbarId (samples.py)
      signal_samples      - samples with signal gen matching (Zbb/Hbb flags, matchedGen_*)
      trigger_sfs         - {era: {"ele": file|None, "mu": file|"pog"|None}}: derived trigger SF
                            files (EleTrigSF / MuTrigSF), "pog" = MUO reference SF, None = 1
      spanet_models       - {name: {path, jet_features, fatjet_features, clean_jets}} (may be empty)
      prefilter_btag_mask - drop events failing btag_mask before stage 2 (analysis selection only)
      dctr_hdamp          - {samples, models: {up, down}}: DCTR hdamp weights (dctr.py); None = off
    '''

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.process_labels = self.workflow_options["process_labels"]
        self.signal_samples = self.workflow_options["signal_samples"]
        self.trigger_sfs = self.workflow_options.get("trigger_sfs", {})
        self.spanet_models = self.workflow_options.get("spanet_models", {})
        self.prefilter_btag_mask = self.workflow_options.get("prefilter_btag_mask", True)
        self.dctr_hdamp = self.workflow_options.get("dctr_hdamp")

    def process_extra_before_skim(self):
        # DCTR hdamp weights and their sums over all generated events (dctr.py; off unless configured)
        if self._isMC:
            apply_dctr_hdamp(self.events, self.output, self._sample, self._dataset, self.dctr_hdamp)

    # ------------------------------------------------------------------ objects
    # (MET filters and the golden JSON are standard skim cuts in the config: eventFlags, goldenJson)
    def apply_object_preselection(self, variation):
        # FatJet.mass = GloParT 2-prong regressed mass, in data and MC
        self.events['FatJet'] = set_fatjet_mass(self.events.FatJet)

        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(self.events.Electron, electron_etaSC, "etaSC")

        self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        self.events['SoftMuonGood'] = soft_lep_sel(self.events, "Muon", self.params)
        self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)
        self.events['SoftElectronGood'] = soft_lep_sel(self.events, "Electron", self.params)

        leptons = ak.with_name(ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1),
                               name='PtEtaPhiMCandidate')
        self.events['LeptonGood'] = leptons[ak.argsort(leptons.pt, ascending=False)]
        self.events['JetGood'], self.jetGoodMask = jet_selection_custom(
            self.events, "Jet", self.params, self._year, "LeptonGood", self.params.object_preselection.btag.tagger)
        self.events['FatJetGood'], self.fatJetGoodMask = jet_selection_custom(
            self.events, "FatJet", self.params, self._year, "LeptonGood", "GloParT")
        self.events['bJetGood'] = bjet_sel(self.events.JetGood, self.params, self._year)
        self.events['qJetGood'] = qjet_sel(self.events.JetGood, self.params, self._year)
        self.events['e_softe'] = lep_softlep_combo(self.events.ElectronGood, self.events.SoftElectronGood)
        self.events['mu_softmu'] = lep_softlep_combo(self.events.MuonGood, self.events.SoftMuonGood)

    def count_objects(self, variation):
        self.events['nMuonGood'] = ak.num(self.events.MuonGood)
        self.events['nElectronGood'] = ak.num(self.events.ElectronGood)
        self.events['nSoftMuonGood'] = ak.num(self.events.SoftMuonGood)
        self.events['nSoftElectronGood'] = ak.num(self.events.SoftElectronGood)
        self.events['nJetGood'] = ak.num(self.events.JetGood)
        self.events['nFatJetGood'] = ak.num(self.events.FatJetGood)
        self.events['nbJetGood'] = ak.num(self.events.bJetGood)
        self.events['nLeptonGood'] = self.events['nMuonGood'] + self.events['nElectronGood']

    def define_common_variables_after_presel(self, variation):
        # b-tag working-point flags of the AK4 jets (UParT)
        for wp, val in self.params.btagging.working_point[self._year]["btagging_WP"].items():
            self.events["JetGood"] = ak.with_field(
                self.events.JetGood, ak.values_astype(self.events.JetGood.btagB > val, int), f"btag_{wp}")
        self.events["LeptonGood"] = ak.with_field(
            self.events.LeptonGood, ak.values_astype(self.events.LeptonGood.pdgId == 11, bool), "is_electron")
        # GloParT XbbVsQCD (btagBB) also as xbbVsQCD, the name used by the Run 2 tools
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, self.events.FatJetGood.btagBB, 'xbbVsQCD')

    # ------------------------------------------------------------------ main hook
    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "btagBB")

        # ---- stage 1
        if self._isMC:
            self.compute_mc_weights()
            match_gen_lep(self.events)
            match_gen_bb_to_zh(self.events)
        self.label_process()
        calc_weight(self.events, self.output, self._dataset, self.params, self._year)

        # ---- prefilter: the expensive part only runs on events of the btag_mask category
        if self.prefilter_btag_mask:
            self.apply_btag_prefilter()

        # ---- stage 2
        fj0 = self.events.FatJetSorted[:, 0]
        self.events['ZH_rap'] = vector.zip({"pt": fj0.pt, "eta": fj0.eta, "phi": fj0.phi, "mass": fj0.mass}).rapidity
        zh_helper(self.events)
        get_controlvars(self.events)
        if self.spanet_models:
            self.onnx_inference()

    # ------------------------------------------------------------------ stage 1 helpers
    def compute_mc_weights(self):
        ev = self.events
        ev['ele_reco_sf'], ev['ele_reco_sfup'], ev['ele_reco_sfdown'] = sf_ele_reco(self.params, ev, self._year)
        ev['ele_id_sf'], ev['ele_id_sfup'], ev['ele_id_sfdown'] = sf_ele_id(self.params, ev, self._year)
        ev['mu_id_sf'], ev['mu_id_sfup'], ev['mu_id_sfdown'], _ = sf_mu(self.params, ev, self._year, 'id')
        ev['mu_iso_sf'], ev['mu_iso_sfup'], ev['mu_iso_sfdown'], _ = sf_mu(self.params, ev, self._year, 'iso')

        # Trigger SFs: our derived files (TrigEff.ipynb), the MUO reference ("pog"), or 1
        trig = self.trigger_sfs.get(self._year, {})
        ele_trig, mu_trig = trig.get("ele"), trig.get("mu")
        if ele_trig:
            ev['ele_trig_sf'], ev['ele_trig_sfup'], ev['ele_trig_sfdown'] = sf_ele_trig(ev, ele_trig)
        else:
            ev['ele_trig_sf'] = ev['ele_trig_sfup'] = ev['ele_trig_sfdown'] = _ones(ev)
        if mu_trig == "pog":
            ev['mu_trig_sf'], ev['mu_trig_sfup'], ev['mu_trig_sfdown'], _ = sf_mu(self.params, ev, self._year, 'trigger')
        elif mu_trig:
            ev['mu_trig_sf'], ev['mu_trig_sfup'], ev['mu_trig_sfdown'] = sf_mu_trig(ev, mu_trig)
        else:
            ev['mu_trig_sf'] = ev['mu_trig_sfup'] = ev['mu_trig_sfdown'] = _ones(ev)

        # b tagging: fixed-WP SFs (UParT, method 1a) need the b-tag efficiency maps; 1 until then
        ev['btag_sf'] = ev['btag_sfup'] = ev['btag_sfdown'] = _ones(ev)

        ev['puWeight'], ev['puWeight_up'], ev['puWeight_down'] = sf_pileup_reweight(self.params, ev, self._year)
        ev['isr'], ev['isr_up'], ev['isr_down'] = sf_partonshower_isr(ev)
        ev['fsr'], ev['fsr_up'], ev['fsr_down'] = sf_partonshower_fsr(ev)
        set_default_dctr_columns(ev)          # 1 where no DCTR weights were computed before the skim
        scale = sf_scaleweights(ev)
        if scale is None:
            scale = (_ones(ev),) * 6
        ev['mu_r_up'], ev['mu_r_down'], ev['mu_f_up'], ev['mu_f_down'], ev['mu_rf_up'], ev['mu_rf_down'] = scale

    def label_process(self):
        """'process' column (+ gen-level flags and top-pT weight where relevant)."""
        if self._sample in self.signal_samples:
            match_gen_sig(self.events, self._sample)
        if ("TTbb" in self._sample) or ("TTTo" in self._sample) or ("ttbb" in self._sample):
            match_gen_tt(self.events, self._sample)          # labels from genTtbarId, top-pT weight
        else:
            if self._sample not in self.process_labels:
                raise ValueError(f"No process label for sample '{self._sample}' (samples.PROCESS_LABELS)")
            self.events['process'] = self.process_labels[self._sample]
            self.events['topptWeight'] = 1
            self.events['topptWeight_Up'] = 1
            self.events['topptWeight_Down'] = 1

    def apply_btag_prefilter(self):
        """Keep only events of the btag_mask category. The category itself still applies the cut,
        so this is purely an optimization (skipped if no event would remain)."""
        self.apply_prefilter(btag_mask)

    def apply_prefilter(self, cut):
        """Keep only events passing `cut` (the category cut, so the output is unchanged)."""
        mask = ak.to_numpy(cut.get_mask(self.events, processor_params=self.params, year=self._year,
                                        sample=self._sample, isMC=self._isMC))
        if mask.all() or not mask.any():
            return
        self.events = self.events[mask]
        self.nEvents_after_presel = self.nevents

    # ------------------------------------------------------------------ SPANet
    def onnx_inference(self):
        models = self.spanet_models
        ev = self.events
        n = len(ev)
        fj0 = ev.FatJetSorted[:, 0]
        lep0 = ev.LeptonGood[:, 0]
        met = ev.PuppiMET

        def jet_features(jets):
            padded = {k: ak.to_numpy(ak.fill_none(ak.pad_none(jets[k], 16, clip=True), 0.))
                      for k in ("pt", "eta", "phi", "btagB")}
            feats = {"pt": np.log(1 + padded["pt"]), "eta": padded["eta"],
                     "sin_phi": np.sin(padded["phi"]), "cos_phi": np.cos(padded["phi"]),
                     "btagB": padded["btagB"]}
            return feats, padded["pt"] != 0

        jet_inputs = {}
        if any(not m.get("clean_jets", False) for m in models.values()):
            jet_inputs[False] = jet_features(ev.JetGood)
        if any(m.get("clean_jets", False) for m in models.values()):
            dr = ev.JetGood.delta_r(ak.firsts(ev.FatJetSorted))
            jet_inputs[True] = jet_features(ev.JetGood[(dr > 0.8) | ak.is_none(dr)])

        lep_data = np.stack([np.log(1 + ak.to_numpy(lep0.pt)), ak.to_numpy(lep0.eta),
                             np.sin(ak.to_numpy(lep0.phi)), np.cos(ak.to_numpy(lep0.phi))],
                            axis=1)[:, None, :].astype(np.float32)
        met_data = np.stack([np.log(1 + ak.to_numpy(met.pt)),
                             np.sin(ak.to_numpy(met.phi)), np.cos(ak.to_numpy(met.phi))],
                            axis=1)[:, None, :].astype(np.float32)
        fatjet_vars = {"pt": np.log(1 + ak.to_numpy(fj0.pt)), "eta": ak.to_numpy(fj0.eta),
                       "sin_phi": np.sin(ak.to_numpy(fj0.phi)), "cos_phi": np.cos(ak.to_numpy(fj0.phi)),
                       "bbtag": ak.to_numpy(fj0.btagBB)}
        mask_global = np.ones((n, 1), dtype=bool)

        for model_name, info in models.items():
            missing = [f for f in info["fatjet_features"] if f not in fatjet_vars]
            if missing:
                raise ValueError(f"SPANet model {model_name}: fatjet features {missing} not available in Run 3")
            jet_vars, jet_mask = jet_inputs[info.get("clean_jets", False)]
            jet_data = np.stack([jet_vars[f] for f in info["jet_features"]], axis=-1).astype(np.float32)
            fatjet_data = np.stack([fatjet_vars[f] for f in info["fatjet_features"]], axis=1)[:, None, :].astype(np.float32)
            inputs = {"Jet_data": jet_data, "Jet_mask": jet_mask,
                      "Met_data": met_data, "Met_mask": mask_global,
                      "Lepton_data": lep_data, "Lepton_mask": mask_global,
                      "FatJet_data": fatjet_data, "FatJet_mask": mask_global}
            if info.get("event_features"):
                # event-level inputs (../common/event_variables.py), on all AK4 jets as in the training
                features = info["event_features"]
                wps = [f[len("n_bjets_"):] for f in features if f.startswith("n_bjets_")]
                values = event_variables(ev.JetGood, ak.firsts(ev.FatJetSorted), ak.firsts(ev.LeptonGood),
                                         wps, info.get("count_clean_jets", False))
                inputs["Event_data"] = np.stack([values[f] for f in features], axis=1)[:, None, :].astype(np.float32)
                inputs["Event_mask"] = mask_global
            probs, encoder = _onnx_session(info["path"]).run(["EVENT/signal", "encoder_event_vector"], inputs)
            if probs.shape[1] != len(SPANET_CLASSES):
                raise ValueError(f"{model_name}: {probs.shape[1]} SPANet classes, expected {SPANET_CLASSES}")
            score = {name: probs[:, i] for name, i in SPANET_CLASSES.items()}
            signal = score["tth"] + score["ttz"]
            self.events[f"spanet_output_{model_name}"] = ak.zip({"signal": signal / probs.sum(axis=1)})
            self.events[f"sig_score_{model_name}"] = signal
            for name, values in score.items():
                self.events[f"{name}_score_{model_name}"] = values
            for i in range(encoder.shape[1]):
                self.events[f"encoderFeature{i + 1}_{model_name}"] = encoder[:, i]
