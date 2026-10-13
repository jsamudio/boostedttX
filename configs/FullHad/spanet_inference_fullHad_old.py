import awkward as ak
import numpy as np
from pocket_coffea.workflows.base import BaseProcessorABC
from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.hist_manager import Axis
from pocket_coffea.lib.objects import (
    jet_correction,
    jet_selection,
    btagging,
)
from object_cleaning_functions import fatjet_sel, bjet_sel, qjet_sel, jet_sel
from custom_cut_functions import sortbyscore
from cand_helperDi import zh_helper
from genmatcher import match_gen_tt, match_gen_sig
from weight_handler import calc_weight, add_weights_to_ttbb
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
        
        # Call base class methods
        zh_helper(self.events)
        
        if self._sample in sig:
            match_gen_sig(self.events, self._sample)
        if ("TTbb" in self._sample) or ("TTTo" in self._sample):
            match_gen_tt(self.events, self._sample)
        if "QCD_HT" in self._sample:
            self.events['process'] = "QCD"
            self.events['topptWeight'] = 1
        
        calc_weight(self.events, self.output, self._dataset, self.params)
        
        if 'TTbb' in self._sample:
            add_weights_to_ttbb(self.events, self._sample)
        
        # Run SPANet inference
        self.onnx_inference(self.workflow_options["spanet_model"])

    def onnx_inference_5feat(self, model_file) -> None:
        """Run SPANet inference - compatible with old 5-feature model.
        
        Jet:    [pt, eta, sin_phi, cos_phi, btag]  → 5 features
        Met:    [pt, eta, sin_phi, cos_phi]         → 4 features
        FatJet: [pt, eta, sin_phi, cos_phi, bbtag]  → 5 features

        Output class indices (EVENT/signal softmax [B, 6]):
            index 0 → unused
            index 1 → tthbb
            index 2 → ttzbb
            index 3 → ttbb
            index 4 → ttlf
            index 5 → qcd
        """

        try:
            worker = get_worker()
        except ValueError:
            worker = None

        if worker is None:
            import onnxruntime as ort
            sess_options = ort.SessionOptions()
            sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            model_session = ort.InferenceSession(
                model_file,
                sess_options=sess_options,
                providers=['CPUExecutionProvider']
            )
        else:
            model_session = worker.data['model_session']

        n_jets_max = 16
        btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]

        # ------------------------------------------------------------------
        # Jet: pt, eta, sin_phi, cos_phi, btag — 5 features, [B, 16, 5]
        # ------------------------------------------------------------------
        pad_dict = {"pt": 0., "phi": 0., "eta": 0., btagging_algorithm: 0.}
        jets_padded = ak.zip({
            key: ak.fill_none(ak.pad_none(self.events.JetGood[key], n_jets_max, clip=True), value)
            for key, value in pad_dict.items()
        })

        jet_pt   = ak.to_numpy(jets_padded.pt)
        jet_mask = (jet_pt != 0)

        jet_data = np.stack([
            np.log(1 + jet_pt),
            ak.to_numpy(jets_padded.eta),
            np.sin(ak.to_numpy(jets_padded.phi)),
            np.cos(ak.to_numpy(jets_padded.phi)),
            ak.to_numpy(jets_padded[btagging_algorithm]),
        ], axis=-1).astype(np.float32)  # [B, 16, 5]

        # ------------------------------------------------------------------
        # MET: pt, eta, sin_phi, cos_phi — 4 features, [B, 1, 4]
        # ------------------------------------------------------------------
        met_eta = (ak.to_numpy(self.events.MET.eta)
                   if hasattr(self.events.MET, 'eta')
                   else np.zeros(len(self.events), dtype=np.float32))

        met_data = np.stack([
            np.log(1 + ak.to_numpy(self.events.MET.pt)),
            met_eta,
            np.sin(ak.to_numpy(self.events.MET.phi)),
            np.cos(ak.to_numpy(self.events.MET.phi)),
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 4]

        met_mask = np.ones((met_data.shape[0], 1), dtype=bool)

        # ------------------------------------------------------------------
        # FatJet: pt, eta, sin_phi, cos_phi, bbtag — 5 features, [B, 1, 5]
        # ------------------------------------------------------------------
        fatjet_data = np.stack([
            np.log(1 + ak.to_numpy(self.events.FatJetSorted[:, 0].pt)),
            ak.to_numpy(self.events.FatJetSorted[:, 0].eta),
            np.sin(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),
            np.cos(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),
            ak.to_numpy(self.events.FatJetSorted[:, 0].xbbVsQCD),
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 5]

        fatjet_mask = np.ones((fatjet_data.shape[0], 1), dtype=bool)

        # ------------------------------------------------------------------
        # Run inference
        # ------------------------------------------------------------------
        outputs = model_session.run(
            input_feed={
                "Jet_data":    jet_data,
                "Jet_mask":    jet_mask,
                "Met_data":    met_data,
                "Met_mask":    met_mask,
                "FatJet_data": fatjet_data,
                "FatJet_mask": fatjet_mask,
            },
            output_names=[
                't1_assignment_probability',
                't2_assignment_probability',
                't1_detection_probability',
                't2_detection_probability',
                'EVENT/signal',
            ]
        )

        t1_assign, t2_assign, t1_detect, t2_detect, event_signal = outputs
        # event_signal: [B, 6]
        # 0=unused, 1=tthbb, 2=ttzbb, 3=ttbb, 4=ttlf, 5=qcd

        self.events["spanet_t1_assignment"] = ak.from_numpy(t1_assign)
        self.events["spanet_t2_assignment"] = ak.from_numpy(t2_assign)
        self.events["spanet_t1_detection"]  = ak.from_numpy(t1_detect)
        self.events["spanet_t2_detection"]  = ak.from_numpy(t2_detect)

        self.events["spanet_output"] = ak.zip({
            "tthbb":  ak.from_numpy(event_signal[:, 1]),
            "ttzbb":  ak.from_numpy(event_signal[:, 2]),
            "ttbb":   ak.from_numpy(event_signal[:, 3]),
            "ttlf":   ak.from_numpy(event_signal[:, 4]),
            "qcd":    ak.from_numpy(event_signal[:, 5]),
            "signal": ak.from_numpy(event_signal[:, 1] + event_signal[:, 2]),
        })

        print(f"SPANet inference complete. "
              f"Mean QCD score: {event_signal[:, 5].mean():.4f}, "
              f"Mean signal score: {(event_signal[:, 1]+event_signal[:, 2]).mean():.4f}")
        
    def onnx_inference(self, model_file) -> None:
        """Run SPANet inference for fully hadronic channel.

        Input features match trainingfile.yaml exactly:
            Jet:    [pt, eta, sin_phi, cos_phi, btag, btag_L, btag_M, btag_H]  → 8 features
            Met:    [pt, eta, sin_phi, cos_phi]                                 → 4 features
            FatJet: [pt, eta, sin_phi, cos_phi, bbtag]                          → 5 features
            Event:  [n_b_inZH]                                                  → 1 feature

        Output class indices (EVENT/signal softmax, shape [B, 6]):
            index 0 → unused (no class 0 in mapping_encoding)
            index 1 → tthbb
            index 2 → ttzbb
            index 3 → ttbb
            index 4 → ttlf
            index 5 → qcd
        """

        try:
            worker = get_worker()
        except ValueError:
            worker = None

        if worker is None:
            import onnxruntime as ort
            sess_options = ort.SessionOptions()
            sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
            model_session = ort.InferenceSession(
                model_file,
                sess_options=sess_options,
                providers=['CPUExecutionProvider']
            )
        else:
            model_session = worker.data['model_session']

        # ------------------------------------------------------------------
        # Verify ONNX input count matches what we expect.
        # Useful for catching model/inference version mismatches early.
        # ------------------------------------------------------------------
        onnx_inputs = {inp.name: inp.shape for inp in model_session.get_inputs()}
        expected_inputs = {'Jet_data', 'Jet_mask', 'Met_data', 'Met_mask',
                           'FatJet_data', 'FatJet_mask', 'Event_data', 'Event_mask'}
        missing = expected_inputs - set(onnx_inputs.keys())
        if missing:
            print(f"WARNING: ONNX missing inputs: {missing}. "
                  f"Export may have been done before Event collection was added. "
                  f"Re-export with: python -m spanet.export <log_dir> spanet.onnx")

        n_jets_max = 16
        btagging_algorithm = self.params.btagging.working_point[self._year]["btagging_algorithm"]

        # ------------------------------------------------------------------
        # Jet features: pt, eta, sin_phi, cos_phi, btag, btag_L, btag_M, btag_H
        # Pad to n_jets_max, clip if more jets exist.
        # ------------------------------------------------------------------
        pad_dict = {"pt": 0., "eta": 0., "phi": 0., btagging_algorithm: 0.,
                    "btag_L": 0., "btag_M": 0., "btag_H": 0.}

        jets_padded = ak.zip({
            key: ak.fill_none(ak.pad_none(self.events.JetGood[key], n_jets_max, clip=True), value)
            for key, value in pad_dict.items()
        })

        jet_pt   = ak.to_numpy(jets_padded.pt)
        jet_mask = (jet_pt != 0)

        jet_data = np.stack([
            np.log(1 + jet_pt),
            ak.to_numpy(jets_padded.eta),
            np.sin(ak.to_numpy(jets_padded.phi)),
            np.cos(ak.to_numpy(jets_padded.phi)),
            ak.to_numpy(jets_padded[btagging_algorithm]),
            ak.to_numpy(jets_padded.btag_L),
            ak.to_numpy(jets_padded.btag_M),
            ak.to_numpy(jets_padded.btag_H),
        ], axis=-1).astype(np.float32)  # [B, 16, 8]

        # ------------------------------------------------------------------
        # MET features: pt, eta, sin_phi, cos_phi  (no phi raw — use sin/cos)
        # Shape: [B, 1, 4]
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
        # FatJet features: pt, eta, sin_phi, cos_phi, bbtag
        # Use leading fat jet (index 0 after sorting by xbbVsQCD).
        # Shape: [B, 1, 5]
        # ------------------------------------------------------------------
        fatjet_data = np.stack([
            np.log(1 + ak.to_numpy(self.events.FatJetSorted[:, 0].pt)),
            ak.to_numpy(self.events.FatJetSorted[:, 0].eta),
            np.sin(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),
            np.cos(ak.to_numpy(self.events.FatJetSorted[:, 0].phi)),
            ak.to_numpy(self.events.FatJetSorted[:, 0].xbbVsQCD),
        ], axis=-1)[:, None, :].astype(np.float32)  # [B, 1, 5]

        fatjet_mask = np.ones((fatjet_data.shape[0], 1), dtype=bool)

        # ------------------------------------------------------------------
        # Event features: n_b_inZH
        # Shape: [B, 1, 1]
        # ------------------------------------------------------------------
        event_data = ak.to_numpy(self.events.n_b_inZH)[:, None, None].astype(np.float32)
        event_mask = np.ones((event_data.shape[0], 1), dtype=bool)

        # ------------------------------------------------------------------
        # Build input feed — only include Event_data if ONNX supports it
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
        else:
            print("WARNING: Event_data not in ONNX — n_b_inZH will not be used. "
                  "Re-export model to include this feature.")

        # ------------------------------------------------------------------
        # Run inference
        # ------------------------------------------------------------------
        output_names = [
            't1_assignment_probability',
            't2_assignment_probability',
            'h_assignment_probability',
            't1_detection_probability',
            't2_detection_probability',
            'h_detection_probability',
            'EVENT/signal',
        ]

        outputs = model_session.run(
            input_feed=input_feed,
            output_names=output_names,
        )

        t1_assign, t2_assign, h_assign, t1_detect, t2_detect, h_detect, event_signal = outputs
        # event_signal shape: [B, 6]
        # index 0 → unused, 1 → tthbb, 2 → ttzbb, 3 → ttbb, 4 → ttlf, 5 → qcd

        self.events["spanet_t1_assignment"] = ak.from_numpy(t1_assign)
        self.events["spanet_t2_assignment"] = ak.from_numpy(t2_assign)
        self.events["spanet_h_assignment"]  = ak.from_numpy(h_assign)
        self.events["spanet_t1_detection"]  = ak.from_numpy(t1_detect)
        self.events["spanet_h_detection"]   = ak.from_numpy(h_detect)
        self.events["spanet_t2_detection"]  = ak.from_numpy(t2_detect)

        self.events["spanet_output"] = ak.zip({
            "tthbb":  ak.from_numpy(event_signal[:, 1]),
            "ttzbb":  ak.from_numpy(event_signal[:, 2]),
            "ttbb":   ak.from_numpy(event_signal[:, 3]),
            "ttlf":   ak.from_numpy(event_signal[:, 4]),
            "qcd":    ak.from_numpy(event_signal[:, 5]),
            # Combined signal score: tthbb + ttzbb probability mass
            "signal": ak.from_numpy(event_signal[:, 1] + event_signal[:, 2]),
        })

        print(f"SPANet inference complete. "
              f"Mean QCD score: {event_signal[:, 5].mean():.4f}, "
              f"Mean signal score: {(event_signal[:, 1] + event_signal[:, 2]).mean():.4f}")
