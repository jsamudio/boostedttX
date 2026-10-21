import pickle

import awkward as ak
import numpy as np
import vector

from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.correction_cache import load_correction_set
from object_cleaning_functions import soft_lep_sel, lep_sel, lep_selTrig_Mu, bjet_sel, qjet_sel, lep_softlep_combo, jet_selection_custom
from custom_cut_functions import sortbyscore, btag_mask
from cand_helper import zh_helper
from genmatcher import match_gen_lep, match_gen_tt, match_gen_sig, match_gen_bb_to_zh
from weight_handler import calc_weight
from pocket_coffea.lib.scale_factors import sf_ele_reco, sf_ele_id, sf_btag, sf_pileup_reweight, sf_L1prefiring
from custom_scale_factors import sf_jet_puId_era, set_fatjet_mass, sf_mu, sf_partonshower_isr, sf_partonshower_fsr, sf_scaleweights, pdfweight_columns, apply_jms_jmr, sf_ele_trig, sf_mu_trig, chunk_rng, btag_jes_systematic, sf_btag_jes
from get_controlvars import get_controlvars
from btag_norm import calib_sums, OTHER_WEIGHTS as BTAG_CALIB_WEIGHTS, VARIATIONS as BTAG_CALIB_VARIATIONS
from theory_norms import fill_theory_sums
from spanet_models import SPANET_CLASSES
from event_variables import event_variables

sig = ['ttHTobb', 'ttHToNonbb', 'TTZToBB', 'TTZToQQ', 'TTZToLLNuNu', 'ttHSMEFT', 'ttZSMEFT']

# Map PocketCoffea year strings to the correctionlib era keys of the PNet mass JSON
PNET_MASS_ERA = {"2016_PreVFP": "2016APV", "2016_PostVFP": "2016", "2017": "2017", "2018": "2018"}

# Legacy process labels -> keys of the Xbb efficiency maps.
# 'old_ttZbb' events are ordinary ttZ real-bb candidates (TTZToQQ events with Z->bb, removed
# downstream in favour of TTZToBB); their tagging efficiency is the ttZ one.
XBB_PROCESS_ALIASES = {'old_ttZbb': 'ttZ', 'ttX': 'TTX'}

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


def _cached(key, loader):
    cache = _worker_cache()
    if key not in cache:
        cache[key] = loader()
    return cache[key]


def _onnx_session(path):
    def load():
        import onnxruntime as ort
        opts = ort.SessionOptions()
        opts.intra_op_num_threads = 1
        opts.inter_op_num_threads = 1
        opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        return ort.InferenceSession(path, sess_options=opts, providers=['CPUExecutionProvider'])
    return _cached(("onnx", path), load)


class ZHbbSpanetProcessor(BaseProcessorABC):
    '''
    Inference workflow. Per chunk and shape variation, after the preselection:

      stage 1 (all preselected events, i.e. without the b-tag requirement)
          object/event scale factors, process labels, Xbb SFs and the sums for the
          b-tag SF normalization (output['btag_calib'], read by btag_norm.py)
      prefilter to the btag_mask category (workflow option 'prefilter_btag_mask')
      stage 2 (btag_mask events only)
          ZH candidate variables, control variables, SPANet inference

    Workflow options:
      pnet_mass_json      - JMS/JMR correctionlib file
      sf_dir              - directory with the Xbb SF, efficiency-map and trigger SF files
      spanet_models       - {name: {path, jet_features, fatjet_features, clean_jets}} (spanet_models.py)
      xbb_wp_cuts         - {era: {L, M, T}} XbbVsQCD thresholds (analysis_years.py)
      btag_calib          - accumulate the b-tag normalization sums (default True; binning in btag_norm.py)
      prefilter_btag_mask - drop events failing btag_mask before stage 2 (default True)
    '''

    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        self.pnet_mass_json = self.workflow_options["pnet_mass_json"]
        self.sf_dir = self.workflow_options["sf_dir"]
        self.xbb_wp = self.workflow_options["xbb_wp_cuts"]
        self.btag_calib = self.workflow_options.get("btag_calib", True)
        self.prefilter_btag_mask = self.workflow_options.get("prefilter_btag_mask", True)
        # Trigger efficiencies (selection "trigeff" in inference_common.py):
        #   trigger_reference_muon - muons = MUO trigger-SF denominator selection (lep_selTrig_Mu:
        #                            tight ID + PFIsoTight, cuts in params/skim_params.yaml) instead of
        #                            the analysis muon (lep_sel); used when the muon trigger is the reference
        #   hlt_columns            - HLT paths exported as columns; paths absent in an era are added as False
        self.trigger_reference_muon = self.workflow_options.get("trigger_reference_muon", False)
        self.hlt_columns = self.workflow_options.get("hlt_columns", [])
        if "spanet_models" not in self.workflow_options:
            raise ValueError("Key `spanet_models` not found in workflow options. Please specify a dictionary of ONNX models.")

    def process_extra_before_skim(self):
        # Theory-weight sums over ALL generated events, for the process norms (theory_norms.py,
        # make_process_norms.py). Runs before the skim, so no event selection is applied.
        if self._isMC:
            fill_theory_sums(self.events, self.output, self._sample, self._dataset)

    # ------------------------------------------------------------------ objects
    # (MET filters and the golden JSON are standard skim cuts in the config: eventFlags, goldenJson)
    def apply_object_preselection(self, variation):
        if self._isMC:
            syst_jms = {"JMSUp": "up", "JMSDown": "down"}.get(variation, "nominal")
            syst_jmr = {"JMRUp": "up", "JMRDown": "down"}.get(variation, "nominal")
            pnet_era = PNET_MASS_ERA.get(self._year, self._year)
            pnet_mass_cset = load_correction_set(self.pnet_mass_json)
            jms_sf = pnet_mass_cset["FatJet_PNet_JMS"].evaluate(pnet_era, syst_jms)
            jmr_sf = pnet_mass_cset["FatJet_PNet_JMR"].evaluate(pnet_era, syst_jmr)
            # The JMAR team used a 10% baseline for the MC resolution in their derivations
            sigma_mc = 0.10
            # Corrections applied to the raw FatJet collection before cuts. The stochastic
            # smearing is seeded per chunk, so every shape variation sees the same random numbers.
            corrected_mass = apply_jms_jmr(self.events, self.events.FatJet, jms_sf, jmr_sf, sigma_mc,
                                           rng=chunk_rng(self.events))
        else:
            # Data: the regressed ParticleNet mass without JMS/JMR (MC-only corrections)
            corrected_mass = self.events.FatJet.particleNet_mass
        self.events['FatJet'] = set_fatjet_mass(self.events.FatJet, corrected_mass)

        # Just needed for plotting
        electron_etaSC = self.events.Electron.eta + self.events.Electron.deltaEtaSC
        self.events["Electron"] = ak.with_field(self.events.Electron, electron_etaSC, "etaSC")

        if self.trigger_reference_muon:
            self.events['MuonGood'] = lep_selTrig_Mu(self.events, "Muon", self.params)
        else:
            self.events['MuonGood'] = lep_sel(self.events, "Muon", self.params)
        self.events['SoftMuonGood'] = soft_lep_sel(self.events, "Muon", self.params)

        # HLT paths of other eras (exported for the trigger efficiencies) are added as False
        for path in self.hlt_columns:
            if path not in self.events.HLT.fields:
                self.events = ak.with_field(self.events, ak.zeros_like(self.events.run, dtype=bool),
                                            where=["HLT", path])
        self.events['ElectronGood'] = lep_sel(self.events, "Electron", self.params)
        self.events['SoftElectronGood'] = soft_lep_sel(self.events, "Electron", self.params)

        leptons = ak.with_name(
                ak.concatenate((self.events.MuonGood, self.events.ElectronGood), axis=1),
                name='PtEtaPhiMCandidate')
        self.events['LeptonGood'] = leptons[ak.argsort(leptons.pt, ascending=False)]
        self.events['JetGood'], self.jetGoodMask = jet_selection_custom(self.events, "Jet", self.params, self._year, "LeptonGood", self.params.object_preselection.btag.tagger)
        self.events['FatJetGood'], self.fatJetGoodMask = jet_selection_custom(self.events, "FatJet", self.params, self._year, "LeptonGood", "PNetMD")
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
        self.events['nFatJetGood2'] = ak.num(self.events.FatJetGood)
        self.events['nbJetGood'] = ak.num(self.events.bJetGood)
        self.events['nLeptonGood'] = (self.events['nMuonGood'] + self.events['nElectronGood'])

    def define_common_variables_after_presel(self, variation):
        super().define_common_variables_before_presel(variation=variation)

        # Labels for btagged jets at the different working points
        for wp, val in self.params.btagging.working_point[self._year]["btagging_WP"].items():
            self.events["JetGood"] = ak.with_field(
                self.events.JetGood,
                ak.values_astype(self.events.JetGood.btagDeepFlavB > val, int),
                f"btag_{wp}"
            )
        self.events["LeptonGood"] = ak.with_field(
            self.events.LeptonGood,
            ak.values_astype(self.events.LeptonGood.pdgId == 11, bool),
            "is_electron"
        )
        fj = self.events.FatJetGood
        xbbVsQCD = fj.particleNetMD_Xbb / (fj.particleNetMD_Xbb + fj.particleNetMD_QCD)
        fj = ak.with_field(fj, xbbVsQCD, 'xbbVsQCD')
        # XbbVsQCD working-point flags (SPANet inputs), per era as in the training workflow
        for wp, cut in self.xbb_wp[self._year].items():
            fj = ak.with_field(fj, ak.values_astype(xbbVsQCD > cut, int), f'xbb_wp_{wp}')
        self.events['FatJetGood'] = fj

    # ------------------------------------------------------------------ main hook
    def process_extra_after_presel(self, variation):
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")

        # ---- stage 1: everything the b-tag normalization needs, on all preselected events
        if self._isMC:
            self.compute_mc_weights(variation)
            match_gen_lep(self.events)
            match_gen_bb_to_zh(self.events)
        self.label_process()
        calc_weight(self.events, self.output, self._dataset, self.params, self._year)
        if self._isMC:
            self.events['xbb_sf'], self.events['xbb_sfup'], self.events['xbb_sfdown'] = self.calculate_xbb_sfs()
            if self.btag_calib:
                self.fill_btag_calib(variation)

        # ---- prefilter: the expensive part only runs on events of the btag_mask category
        if self.prefilter_btag_mask:
            self.apply_btag_prefilter()

        # ---- stage 2
        fj0 = self.events.FatJetSorted[:, 0]
        self.events['ZH_rap'] = vector.zip({"pt": fj0.pt, "eta": fj0.eta, "phi": fj0.phi, "mass": fj0.mass}).rapidity
        zh_helper(self.events)
        get_controlvars(self.events)
        self.onnx_inference()

    # ------------------------------------------------------------------ stage 1 helpers
    def compute_mc_weights(self, variation):
        ev = self.events
        ev['ele_reco_sf'], ev['ele_reco_sfup'], ev['ele_reco_sfdown'] = sf_ele_reco(self.params, ev, self._year)
        ev['ele_id_sf'], ev['ele_id_sfup'], ev['ele_id_sfdown'] = sf_ele_id(self.params, ev, self._year)
        ev['ele_trig_sf'], ev['ele_trig_sfup'], ev['ele_trig_sfdown'] = sf_ele_trig(ev, f"{self.sf_dir}/ele_trig_sf_{self._year}.json.gz")
        # muon reco / ID / iso SFs: total up/down, and the stat and syst parts (mu_<key>_sf<part>)
        for key in ['id', 'iso', 'reco']:
            ev[f'mu_{key}_sf'], ev[f'mu_{key}_sfup'], ev[f'mu_{key}_sfdown'], parts = sf_mu(self.params, ev, self._year, key)
            for part, values in parts.items():
                ev[f'mu_{key}_sf{part}'] = values
        if self.trigger_reference_muon:
            # Reference muon (trigger efficiencies, pass 1): central MUO trigger SF, defined for
            # tight ID + PFIsoTight muons (key 'trigger' in params/lepton_scale_factors.yaml)
            ev['mu_trig_sf'], ev['mu_trig_sfup'], ev['mu_trig_sfdown'], _ = sf_mu(self.params, ev, self._year, 'trigger')
        else:
            ev['mu_trig_sf'], ev['mu_trig_sfup'], ev['mu_trig_sfdown'] = sf_mu_trig(ev, f"{self.sf_dir}/muon_trig_sf_{self._year}.json.gz")

        # b-tag shape SFs: btag_sf, btag_sf<src>, btag_sf<src>_up, btag_sf<src>_down
        btag_variations = self.params.systematic_variations.weight_variations.sf_btag[self._year]
        btag_sf = sf_btag(self.params, ev.JetGood, self._year, njets=ev.nJetGood,
                          variations=['central'] + btag_variations)
        ev['btag_sf'] = btag_sf['central'][0]
        for src in btag_variations:
            ev[f'btag_sf{src}'], ev[f'btag_sf{src}_up'], ev[f'btag_sf{src}_down'] = btag_sf[src][:3]
        # AK4 JES shape variations: the central b-tag SF is replaced by the JES-varied one
        jes_syst = btag_jes_systematic(variation)
        if jes_syst is not None:
            ev['btag_sf'] = sf_btag_jes(self.params, ev.JetGood, self._year, ev.nJetGood, jes_syst)

        ev['puWeight'], ev['puWeight_up'], ev['puWeight_down'] = sf_pileup_reweight(self.params, ev, self._year)
        # Pileup jet ID SF (JME jmar.json, PUJetID_eff): AK4 jets with pT < maxpt that pass the PU ID WP
        # and are matched to a GenJet; WP / pT limit from object_preselection.yaml (Jet: puId), PU ID bit
        # value of the WP per era (custom_scale_factors.sf_jet_puId_era)
        ev['jet_puId_sf'], ev['jet_puId_sfup'], ev['jet_puId_sfdown'] = sf_jet_puId_era(self.params, ev.JetGood, self._year, ev.nJetGood)
        ev['isr'], ev['isr_up'], ev['isr_down'] = sf_partonshower_isr(ev)
        ev['fsr'], ev['fsr_up'], ev['fsr_down'] = sf_partonshower_fsr(ev)
        ev['mu_r_up'], ev['mu_r_down'], ev['mu_f_up'], ev['mu_f_down'], ev['mu_rf_up'], ev['mu_rf_down'] = sf_scaleweights(ev)
        # PDF / alpha_s members pdfweight_0..102 (exported for signal and tt+jets, outvars.pdf_vars)
        for name, values in pdfweight_columns(ev).items():
            ev[name] = values
        ev['l1prefire_sf'], ev['l1prefire_sfup'], ev['l1prefire_sfdown'] = sf_L1prefiring(ev)

    def label_process(self):
        """'process' column (+ gen-level flags and top-pT weight where relevant)."""
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        elif ("TTbb" in self._sample) or ("TTTo" in self._sample) or ("ttbb" in self._sample):
            match_gen_tt(self.events, self._sample)
        else:
            self.events['process'] = 'VJets' if "JetsToL" in self._sample else self._sample
            self.events['topptWeight'] = 1
            self.events['topptWeight_Up'] = 1
            self.events['topptWeight_Down'] = 1

    def calculate_xbb_sfs(self):
        """Per-event Xbb SF (nominal, up, down): multi-WP method for real (2 gen b) candidates,
        data-driven fake SFs otherwise. Returns numpy arrays."""
        era = self._year
        real_eval = load_correction_set(f"{self.sf_dir}/SFBDT_Xbb_{era}.json")[f"SFBDT_Xbb_SF_{era}"]
        fake_eval = load_correction_set(f"{self.sf_dir}/Fake_SF_Xbb_{era}.json")[f"Fake_SF_Xbb_{era}"]

        def load_effmaps():
            with open(f"{self.sf_dir}/efficiency_maps/efficiency_maps_{era}.pkl", 'rb') as f:
                return pickle.load(f)
        eff_maps = _cached(("xbb_eff", era), load_effmaps)

        fatjets = self.events.FatJetSorted[:, 0]
        pt = ak.to_numpy(fatjets.pt).astype(np.float64)
        msd = ak.to_numpy(fatjets.mass).astype(np.float64)
        score = ak.to_numpy(fatjets.xbbVsQCD)
        genmatch = ak.to_numpy(self.events.xbb_genmatch)
        n_events = len(pt)
        if n_events == 0:
            empty = np.ones(0)
            return empty, empty, empty

        # Mutually exclusive tagging categories
        thr = eff_maps['thresholds']
        is_hp = score >= thr['HP']
        is_mp = (score >= thr['MP']) & (score < thr['HP'])
        is_lp = (score >= thr['LP']) & (score < thr['MP'])
        is_fail = score < thr['LP']

        pt_idx = np.clip(np.digitize(pt, eff_maps['pt_bins']) - 1, 0, len(eff_maps['pt_bins']) - 2)
        msd_idx = np.clip(np.digitize(msd, eff_maps['msd_bins']) - 1, 0, len(eff_maps['msd_bins']) - 2)

        # Efficiencies per event from the event's own process label (labels can differ
        # within a sample, e.g. TTZToQQ -> 'ttZ' / 'old_ttZbb')
        process_arr = ak.to_numpy(self.events['process']).copy()
        for old_label, new_label in XBB_PROCESS_ALIASES.items():
            process_arr[process_arr == old_label] = new_label
        per_process = eff_maps['real']['per_process']
        unknown_mask = ~np.isin(process_arr, list(per_process.keys()))
        if np.any(unknown_mask):
            print(f"WARNING: {int(unknown_mask.sum())} event(s) with process label(s) "
                  f"{list(np.unique(process_arr[unknown_mask]))} not in the Xbb efficiency maps "
                  f"(known: {sorted(per_process)}); xbb_sf = 1 for these events.")

        eff = {wp: np.ones(n_events) for wp in ("HP", "MP", "LP")}
        for proc, maps in per_process.items():
            m = process_arr == proc
            if np.any(m):
                for wp in eff:
                    eff[wp][m] = maps[wp]['eff'][pt_idx[m], msd_idx[m]]

        out = []
        for syst in ("nominal", "up", "down"):
            # Real candidates (method 1a, multi-WP): exclusive-bin weights from the cumulative SFs
            sf_hp = real_eval.evaluate("HP", syst, pt)
            sf_mp = real_eval.evaluate("MP", syst, pt)
            sf_lp = real_eval.evaluate("LP", syst, pt)
            num_mp = sf_mp * eff["MP"] - sf_hp * eff["HP"]
            den_mp = eff["MP"] - eff["HP"]
            w_mp = np.where(num_mp < 0, 1.0, num_mp / np.where(den_mp > 0, den_mp, 1.0))
            num_lp = sf_lp * eff["LP"] - sf_mp * eff["MP"]
            den_lp = eff["LP"] - eff["MP"]
            w_lp = np.where(num_lp < 0, 1.0, num_lp / np.where(den_lp > 0, den_lp, 1.0))
            real_w = np.select([is_hp, is_mp, is_lp], [sf_hp, w_mp, w_lp], default=1.0)

            # Fake candidates (data driven), per tagging category
            fake_w = np.ones(n_events)
            for wp, sel in (("HP", is_hp), ("MP", is_mp), ("LP", is_lp), ("Fail", is_fail)):
                if np.any(sel):
                    fake_w[sel] = fake_eval.evaluate(wp, syst, pt[sel], msd[sel])

            final = np.where(genmatch, real_w, fake_w)
            final[unknown_mask] = 1.0
            out.append(final)
        return tuple(out)

    def fill_btag_calib(self, variation):
        """Accumulate the per-jet-bin sums for btag_norm.py in output['btag_calib'][dataset][sample][variation].
        Sums are without the per-dataset normalization (applied in btag_norm.py)."""
        ev = self.events
        w0 = np.sign(ak.to_numpy(ev.genWeight))
        for col in BTAG_CALIB_WEIGHTS:
            if col in ev.fields:
                w0 = w0 * ak.to_numpy(ev[col])
        sfs = {var: ak.to_numpy(ev[col]) for var, col in BTAG_CALIB_VARIATIONS.items() if col in ev.fields}
        njets = ak.to_numpy(ev.nJetGood)

        if self._hasSubsamples:
            subs = self._subsamples[self._sample]
            subs.prepare(events=ev, processor_params=self.params, year=self._year,
                         sample=self._sample, isMC=self._isMC)
            masks = [(f"{self._sample}__{name}", ak.to_numpy(mask)) for name, mask in subs.get_masks()]
        else:
            masks = [(self._sample, np.ones(len(ev), dtype=bool))]

        out = self.output.setdefault("btag_calib", {}).setdefault(self._dataset, {})
        for key, m in masks:
            out.setdefault(key, {})[variation] = calib_sums(njets[m], w0[m], {v: s[m] for v, s in sfs.items()})

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
        models = self.workflow_options["spanet_models"]
        if not models:          # e.g. the trigger-efficiency groups
            return
        ev = self.events
        n = len(ev)
        fj0 = ev.FatJetSorted[:, 0]
        lep0 = ev.LeptonGood[:, 0]

        def jet_features(jets):
            padded = {k: ak.to_numpy(ak.fill_none(ak.pad_none(jets[k], 16, clip=True), 0.))
                      for k in ("pt", "eta", "phi", "btagB")}
            feats = {
                "pt": np.log(1 + padded["pt"]),
                "eta": padded["eta"],
                "sin_phi": np.sin(padded["phi"]),
                "cos_phi": np.cos(padded["phi"]),
                "btagB": padded["btagB"],
            }
            return feats, padded["pt"] != 0

        jet_inputs = {}
        if any(not m.get("clean_jets", False) for m in models.values()):
            jet_inputs[False] = jet_features(ev.JetGood)
        if any(m.get("clean_jets", False) for m in models.values()):
            # AK4 jets outside the cone of the ZH candidate (as in the training workflow)
            dr = ev.JetGood.delta_r(ak.firsts(ev.FatJetSorted))
            jet_inputs[True] = jet_features(ev.JetGood[(dr > 0.8) | ak.is_none(dr)])

        lep_data = np.stack([np.log(1 + ak.to_numpy(lep0.pt)), ak.to_numpy(lep0.eta),
                             np.sin(ak.to_numpy(lep0.phi)), np.cos(ak.to_numpy(lep0.phi))],
                            axis=1)[:, None, :].astype(np.float32)
        met_data = np.stack([np.log(1 + ak.to_numpy(ev.MET.pt)),
                             np.sin(ak.to_numpy(ev.MET.phi)), np.cos(ak.to_numpy(ev.MET.phi))],
                            axis=1)[:, None, :].astype(np.float32)
        fatjet_vars = {
            "pt": np.log(1 + ak.to_numpy(fj0.pt)),
            "eta": ak.to_numpy(fj0.eta),
            "sin_phi": np.sin(ak.to_numpy(fj0.phi)),
            "cos_phi": np.cos(ak.to_numpy(fj0.phi)),
            "bbtag": ak.to_numpy(fj0.xbbVsQCD),
            "xbb_wp_L": ak.to_numpy(fj0.xbb_wp_L),
            "xbb_wp_M": ak.to_numpy(fj0.xbb_wp_M),
            "xbb_wp_T": ak.to_numpy(fj0.xbb_wp_T),
        }
        mask_global = np.ones((n, 1), dtype=bool)

        for model_name, info in models.items():
            jet_vars, jet_mask = jet_inputs[info.get("clean_jets", False)]
            jet_data = np.stack([jet_vars[f] for f in info["jet_features"]], axis=-1).astype(np.float32)
            fatjet_data = np.stack([fatjet_vars[f] for f in info["fatjet_features"]], axis=1)[:, None, :].astype(np.float32)
            inputs = {
                "Jet_data": jet_data, "Jet_mask": jet_mask,
                "Met_data": met_data, "Met_mask": mask_global,
                "Lepton_data": lep_data, "Lepton_mask": mask_global,
                "FatJet_data": fatjet_data, "FatJet_mask": mask_global,
            }
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
