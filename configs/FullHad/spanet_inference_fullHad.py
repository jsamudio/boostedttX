import awkward as ak
import numpy as np
from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_selection,
    btagging,
)
from object_cleaning_functions import fatjet_sel, bjet_sel, qjet_sel, jet_sel
from custom_cut_functions import sortbyscore
from cand_helperDi import zh_helper
from genmatcher import match_gen_tt, match_gen_sig
from weight_handler import calc_weight
#,add_weights_to_ttbb
from coffea.analysis_tools import PackedSelection
from pocket_coffea.lib.parton_provenance import *
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.scale_factors import sf_btag
from dask.distributed import get_worker

# Import the base class
from zhbbworkflow_spanet import ZHbbBaseProcessor

sig = ['ttHTobb', 'ttHToNonbb','TTZToBB', 'TTZToQQ', 'TTZToLLNuNu']

class ZHbbSpanetProcessor(ZHbbBaseProcessor):
    def __init__(self, cfg: Configurator):
        super().__init__(cfg)
        self.dr_min = self.workflow_options["parton_jet_min_dR"]
        self.dr_min_postfsr = self.workflow_options.get("parton_jet_min_dR_postfsr", 1.)
        if not "spanet_model" in self.workflow_options:
            raise ValueError("Key `spanet_model` not found in workflow options. Please specify the path to the ONNX model.")
        elif not self.workflow_options["spanet_model"].endswith(".onnx"):
            raise ValueError("Key `spanet_model` should be the path of an ONNX model.")
        
    def apply_object_preselection(self, variation):
        """Apply object preselection for fully hadronic channel - no leptons!"""
        
        # Jet selection - no lepton cleaning needed for FH
        self.events['JetGood'], self.jetGoodMask = jet_sel(self.events, "Jet", self.params, "")
        
        # Fat jet selection
        self.events['FatJetGood'] = fatjet_sel(self.events, self.params, "")
        
        # b-jet and light-jet categorization
        self.events['bJetGood'] = bjet_sel(self.events.JetGood, self.params)
        self.events['qJetGood'] = qjet_sel(self.events.JetGood, self.params)
        
        # Add columns for SPANet
        self.events['JetGood_pt'] = self.events.JetGood.pt
        self.events['JetGood_eta'] = self.events.JetGood.eta
        self.events['JetGood_phi'] = self.events.JetGood.phi

    def define_common_variables_after_presel(self, variation):
        # Call base class method first
        super().define_common_variables_after_presel(variation)
        
        # Compute Xbb/(Xbb+QCD) score for fat jets
        xbbVsQCD = self.events.FatJetGood.particleNetMD_Xbb / (
            self.events.FatJetGood.particleNetMD_Xbb + self.events.FatJetGood.particleNetMD_QCD
        )
        self.events['FatJetGood'] = ak.with_field(self.events.FatJetGood, xbbVsQCD, 'xbbVsQCD')
        
        # Sort fat jets by bbtag score
        self.events['FatJetSorted'] = sortbyscore(self.events.FatJetGood, "xbbVsQCD")

    def count_objects(self, variation):
        self.events['nJetGood'] = ak.num(self.events.JetGood)
        self.events['nFatJetGood'] = ak.num(self.events.FatJetGood)
        self.events['nbJetGood'] = ak.num(self.events.bJetGood)
        self.events['nqJetGood'] = ak.num(self.events.qJetGood)
        
        n_events = len(self.events)
        zeros = ak.zeros_like(self.events.event, dtype=np.float32)
        
        # Mimic lep_softlep_combo output: flat per-event record, all zeros
        dummy_dilep = ak.zip({
            "pt":     zeros,
            "eta":    zeros,
            "phi":    zeros,
            "mass":   zeros,
            "charge": zeros,
            "deltaR": ak.full_like(zeros, -1.0),
        }, with_name="PtEtaPhiMCandidate")
        
        self.events['e_softe']  = dummy_dilep
        self.events['mu_softmu'] = dummy_dilep
        
        self.events['nMuonGood']        = ak.zeros_like(self.events.event, dtype=int)
        self.events['nElectronGood']    = ak.zeros_like(self.events.event, dtype=int)
        self.events['nSoftMuonGood']    = ak.zeros_like(self.events.event, dtype=int)
        self.events['nSoftElectronGood'] = ak.zeros_like(self.events.event, dtype=int)
        self.events['nLeptonGood']      = ak.zeros_like(self.events.event, dtype=int)
    
    def process_extra_after_presel(self, variation):
        """Main processing after preselection"""
        self.events['topptWeight'] = ak.ones_like(
            self.events.event, dtype=np.float32
        )        
        # Call base class methods
        zh_helper(self.events)
        
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        if ("TTbb" in self._sample) or ("TTTo" in self._sample):
            match_gen_tt(self.events, self._sample)
        if "QCD_HT" in self._sample:
            self.events['process'] = "QCD"

        
        calc_weight(self.events, self.output, self._dataset, self.params)

        self.events['HT']        = ak.sum(self.events.JetGood.pt, axis=1)
        self.events['nJetGood']  = ak.num(self.events.JetGood)
        self.events['nbJetGood'] = ak.num(self.events.bJetGood)
        # Run SPANet inference
        self.onnx_inference(self.workflow_options["spanet_model"])
    def compute_top_masses(self, t1_assign, t2_assign):
        """
        Compute reconstructed top masses from SPANet assignment probabilities.
        
        t1_assign: [B, N, N, N] — probability for each (b, q1, q2) jet combo for t1
        t2_assign: [B, N, N, N] — same for t2
        
        Returns:
        t1_mass: [B] — reconstructed t1 mass in GeV
        t2_mass: [B] — reconstructed t2 mass in GeV  
        t1_w_mass: [B] — reconstructed W mass from t1 (q1+q2)
        t2_w_mass: [B] — reconstructed W mass from t2 (q1+q2)
        """
        n_jets_max = 16
        B = t1_assign.shape[0]
        
        # ── Pad jet 4-vectors to fixed size ──────────────────────────────
        pad_dict = {
            "pt": 0., "eta": 0., "phi": 0., "mass": 0.
        }
        jets_padded = ak.zip({
            key: ak.fill_none(
                ak.pad_none(self.events.JetGood[key], n_jets_max, clip=True), val
            )
            for key, val in pad_dict.items()
        })
        
        # Convert to Cartesian 4-momenta [B, N]
        pt   = ak.to_numpy(jets_padded.pt).astype(np.float32)    # [B, N]
        eta  = ak.to_numpy(jets_padded.eta).astype(np.float32)
        phi  = ak.to_numpy(jets_padded.phi).astype(np.float32)
        mass = ak.to_numpy(jets_padded.mass).astype(np.float32)
        
        px = pt * np.cos(phi)                                     # [B, N]
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        E  = np.sqrt(pt**2 * np.cosh(eta)**2 + mass**2)
        
        def top_mass_from_assignment(assign):
            """
            assign: [B, N, N, N] — (b_idx, q1_idx, q2_idx) probabilities
            Returns: top_mass [B], w_mass [B]
            """
            # Find best assignment: argmax over flattened N³ combinations
            flat = assign.reshape(B, -1)           # [B, N³]
            best = flat.argmax(axis=-1)            # [B]
            
            b_idx  = best // (n_jets_max * n_jets_max)
            q1_idx = (best % (n_jets_max * n_jets_max)) // n_jets_max
            q2_idx = best % n_jets_max
            
            # Gather 4-momenta for best assignment
            b_range = np.arange(B)
            
            # Top: b + q1 + q2
            top_E  = E[b_range, b_idx]  + E[b_range, q1_idx]  + E[b_range, q2_idx]
            top_px = px[b_range, b_idx] + px[b_range, q1_idx] + px[b_range, q2_idx]
            top_py = py[b_range, b_idx] + py[b_range, q1_idx] + py[b_range, q2_idx]
            top_pz = pz[b_range, b_idx] + pz[b_range, q1_idx] + pz[b_range, q2_idx]
            
            top_mass2 = top_E**2 - top_px**2 - top_py**2 - top_pz**2
            top_mass  = np.sqrt(np.clip(top_mass2, 0, None))
            
            # W: q1 + q2 only
            w_E  = E[b_range, q1_idx]  + E[b_range, q2_idx]
            w_px = px[b_range, q1_idx] + px[b_range, q2_idx]
            w_py = py[b_range, q1_idx] + py[b_range, q2_idx]
            w_pz = pz[b_range, q1_idx] + pz[b_range, q2_idx]
            
            w_mass2 = w_E**2 - w_px**2 - w_py**2 - w_pz**2
            w_mass  = np.sqrt(np.clip(w_mass2, 0, None))
            
            return top_mass, w_mass, top_E, top_px, top_py, top_pz
        
        t1_mass, t1_w_mass, t1_E, t1_px, t1_py, t1_pz = \
            top_mass_from_assignment(t1_assign)
        t2_mass, t2_w_mass, t2_E, t2_px, t2_py, t2_pz = \
            top_mass_from_assignment(t2_assign)
        tt_E  = t1_E  + t2_E
        tt_px = t1_px + t2_px
        tt_py = t1_py + t2_py
        tt_pz = t1_pz + t2_pz
        tt_mass2 = tt_E**2 - tt_px**2 - tt_py**2 - tt_pz**2
        tt_mass  = np.sqrt(np.clip(tt_mass2, 0, None))

        return t1_mass, t2_mass, t1_w_mass, t2_w_mass, tt_mass

    def onnx_inference(self, model_file) -> None:
        """Run SPANet inference for fully hadronic channel - version 11.
        
        Run SPANet inference for fully hadronic channel - version 19.
        
        Jet:    [MASK, pt, eta, sin_phi, cos_phi, btag, btag_L, btag_M, btag_H,
         btagDeepFlavCvL, btagDeepFlavCvB, btagDeepFlavQG,
         pt_frac, dR_fatjet]  → 14 features + MASK = 15
        Met:    [pt, eta, sin_phi, cos_phi]                        → 4 features
        FatJet: [pt, eta, sin_phi, cos_phi, bbtag,
         deepTagMD_ZHbbvsQCD, deepTagMD_TvsQCD,
         particleNet_TvsQCD]                               → 8 features
        Event:  [n_b_inZH, nJetGood, nbJetGood, HT]               → 4 features
        Jet:    [MASK, pt, eta, sin_phi, cos_phi, btag, btag_L, btag_M, btag_H,
             hadronFlavour, pt_frac, dR_fatjet]  → 12 features + MASK = 13
        Met:    [pt, eta, sin_phi, cos_phi]           → 4 features
        FatJet: [pt, eta, sin_phi, cos_phi, bbtag]    → 5 features
        Event:  [n_b_inZH, nJetGood, nbJetGood, HT]  → 4 features
        """
        try:
            worker = get_worker()
        except ValueError:
            worker = None
        
        if worker is None:
            import onnxruntime as ort
            sess_options = ort.SessionOptions()
            sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            sess_options.intra_op_num_threads = 4   
            sess_options.inter_op_num_threads = 2   
            model_session = ort.InferenceSession(
                model_file,
                sess_options=sess_options,
                providers=['CPUExecutionProvider']
            )
        else:
            model_session = worker.data['model_session']

        onnx_inputs = {inp.name: inp.shape for inp in model_session.get_inputs()}

        n_jets_max = 16
        btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]
        
        # ------------------------------------------------------------------
        # Jet features: 13 total (MASK + 12 features)
        # ------------------------------------------------------------------
        pad_dict = {
            "pt": 0., "eta": 0., "phi": 0.,
            btagging_algorithm: 0.,
            "btag_L": 0., "btag_M": 0., "btag_H": 0.,
            "btagDeepFlavCvL": 0.,
            "btagDeepFlavCvB": 0.,
            "btagDeepFlavQG":  0.,
            "pt_frac": 0., "dR_fatjet": 0.
        }        
        jets_padded = ak.zip({
            key: ak.fill_none(
                ak.pad_none(self.events.JetGood[key], n_jets_max, clip=True), value
            )
            for key, value in pad_dict.items()
        })

        jet_pt   = ak.to_numpy(jets_padded.pt)
        jet_mask = (jet_pt != 0)

        jet_data = np.stack([
            jet_mask.astype(np.float32),                           # MASK
            np.log(1 + jet_pt),                                    # pt
            ak.to_numpy(jets_padded.eta),                          # eta
            np.sin(ak.to_numpy(jets_padded.phi)),                  # sin_phi
            np.cos(ak.to_numpy(jets_padded.phi)),                  # cos_phi
            ak.to_numpy(jets_padded[btagging_algorithm]),          # btag
            ak.to_numpy(jets_padded.btag_L),                       # btag_L
            ak.to_numpy(jets_padded.btag_M),                       # btag_M
            ak.to_numpy(jets_padded.btag_H),                       # btag_H
            ak.to_numpy(jets_padded.btagDeepFlavCvL),              # CvL ← NEW
            ak.to_numpy(jets_padded.btagDeepFlavCvB),              # CvB ← NEW
            ak.to_numpy(jets_padded.btagDeepFlavQG),               # QG  ← NEW
            np.log(1 + ak.to_numpy(jets_padded.pt_frac)),          # pt_frac
            ak.to_numpy(jets_padded.dR_fatjet),                    # dR_fatjet
        ], axis=-1).astype(np.float32)  # [B, 16, 14]
        # ------------------------------------------------------------------
        # MET features — unchanged
        # ------------------------------------------------------------------
        met_data = np.stack([
            np.log(1 + ak.to_numpy(self.events.MET.pt)),
            ak.to_numpy(self.events.MET.eta) if hasattr(self.events.MET, 'eta')
            else np.zeros(len(self.events), dtype=np.float32),
            np.sin(ak.to_numpy(self.events.MET.phi)),
            np.cos(ak.to_numpy(self.events.MET.phi)),
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 4]
        
        met_mask = np.ones((met_data.shape[0], 1), dtype=bool)

        # ------------------------------------------------------------------
        # FatJet features — unchanged
        # ------------------------------------------------------------------
        fatjet_data = np.stack([
            np.log(1 + ak.to_numpy(self.events.FatJetSorted[:, 0].pt)),          # pt
            ak.to_numpy(self.events.FatJetSorted[:, 0].eta),                      # eta
            np.sin(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),              # sin_phi
            np.cos(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),              # cos_phi
            ak.to_numpy(self.events.FatJetSorted[:, 0].xbbVsQCD),                 # bbtag
            ak.to_numpy(self.events.FatJetSorted[:, 0].deepTagMD_ZHbbvsQCD),      # ZHbb ← NEW
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 8]        
        fatjet_mask = np.ones((fatjet_data.shape[0], 1), dtype=bool)
        
        # ------------------------------------------------------------------
        # Event features: n_b_inZH, nJetGood, nbJetGood, HT — 4 features
        # ------------------------------------------------------------------
        event_data = np.stack([
            ak.to_numpy(self.events.n_b_inZH).astype(np.float32),
            ak.to_numpy(self.events.nJetGood).astype(np.float32),
            ak.to_numpy(self.events.nbJetGood).astype(np.float32),
            np.log(1 + ak.to_numpy(self.events.HT)),  
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 4]
        
        event_mask = np.ones((event_data.shape[0], 1), dtype=bool)
        
        # ------------------------------------------------------------------
        # Build input feed
        # ------------------------------------------------------------------
        input_feed = {
            "Jet_data":    jet_data,
            "Jet_mask":    jet_mask,
            "Met_data":    met_data,
            "Met_mask":    met_mask,
            "FatJet_data": fatjet_data,
            "FatJet_mask": fatjet_mask,
        }
        
        if 'Event_data' in onnx_inputs:
            input_feed["Event_data"] = event_data
            input_feed["Event_mask"] = event_mask
            
        # ------------------------------------------------------------------
        # Run inference
        # ------------------------------------------------------------------

        outputs = model_session.run(
            input_feed=input_feed,
            output_names=[
                't1_assignment_probability',
                't2_assignment_probability',
                't1_detection_probability',
                't2_detection_probability',
                'EVENT/signal',
            ]
        )
        
        t1_assign, t2_assign, t1_detect, t2_detect, event_signal = outputs

        self.events["spanet_t1_assignment"] = ak.from_numpy(t1_assign)
        self.events["spanet_t2_assignment"] = ak.from_numpy(t2_assign)
        self.events["spanet_t1_detection"]  = ak.from_numpy(t1_detect)
        self.events["spanet_t2_detection"]  = ak.from_numpy(t2_detect)

        t1_assign, t2_assign, t1_detect, t2_detect, event_signal = outputs

        # ── Compute reconstructed top and W masses ──────────────────────────
        t1_mass, t2_mass, t1_w_mass, t2_w_mass, tt_mass = self.compute_top_masses(
            t1_assign, t2_assign
        )
        
        # Store masses as event fields
        self.events["spanet_t1_mass"]   = ak.from_numpy(t1_mass)
        self.events["spanet_t2_mass"]   = ak.from_numpy(t2_mass)
        self.events["spanet_t1_w_mass"] = ak.from_numpy(t1_w_mass)
        self.events["spanet_t2_w_mass"] = ak.from_numpy(t2_w_mass)
        self.events["spanet_tt_mass"]   = ak.from_numpy(tt_mass)
        
        self.events["spanet_output"] = ak.zip({
            "tthbb":  ak.from_numpy(event_signal[:, 0]),
            "ttzbb":  ak.from_numpy(event_signal[:, 1]),
            "ttbb":   ak.from_numpy(event_signal[:, 2]),
            "ttc":    ak.from_numpy(event_signal[:, 3]),
            "ttlf":   ak.from_numpy(event_signal[:, 4]),
            "qcd":    ak.from_numpy(event_signal[:, 5]),
            "signal": ak.from_numpy(event_signal[:, 0] + event_signal[:, 1]),
        })
        
        print(f"SPANet inference complete. "
              f"Mean QCD score: {event_signal[:, 4].mean():.4f}, "
              f"Mean signal score: {(event_signal[:, 0] + event_signal[:, 1]).mean():.4f}")
