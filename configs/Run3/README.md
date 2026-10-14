# Run 3 semileptonic boosted ttZ/ttH analysis (2024)

Same structure as `configs/Run2`, with the Run 3 object definitions. The previous 2024 setup
(`configs/SemiLeptonic_2024`) is kept unchanged for comparison until this one is validated.

## Object definitions

| | 2024 |
|---|---|
| Data / MC | NanoAODv15, Summer24 MC, runs C-I (109.95 fb^-1, golden JSON) |
| MET | PUPPI MET; pocket-coffea's MET calibrator recomputes the type-1 PUPPI MET from `RawPuppiMET` with the recommended JECs (data and MC) |
| JECs | pocket-coffea defaults: `Summer24Prompt24_V5`, JER `Summer24Prompt24_JRV2`, reapplied to data and MC (AK4 PUPPI, CorrT1METJet, AK8) |
| AK4 b tagging | UParTAK4, medium WP from pocket-coffea's btagging parameters (`object_preselection.btag.wp`) |
| AK8 tagger | GloParT XbbVsQCD = Xbb / (Xbb + QCD); the ZH candidate is the leading fatjet in it |
| Fatjet mass | GloParT 2-prong regressed mass (`globalParT3_massCorrX2p` x raw mass), in data and MC, for the mass window, `ZH_M`, the columns and 4-vectors (`custom_scale_factors.set_fatjet_mass`) |
| Triggers | Ele30_WPTight_Gsf, Ele115_CaloIdVT_GsfTrkIdT / IsoMu24, Mu50, HighPtTkMu100, CascadeMu100; per primary dataset in data (Muon PD: muon paths; EGamma PD: e/gamma paths and no muon path) |
| Leptons | electrons MVA noIso WP90, muons medium ID, mini-isolation (`params/object_preselection.yaml`) |

Process labels come from `samples.PROCESS_LABELS` (ttZ = TTZToQQ incl. Z->bb, TTLL, TTNuNu;
ttH = ttHTobb, ttHToNonbb); ttbar samples are split by genTtbarId.

## Not derived yet (bottom-up)

1. **Trigger SFs** (next): `inference_conf_trigeff.py` (ttbar e-mu, no trigger requirement) +
   `TrigEff.ipynb` (orthogonal method, as in Run 2), two passes:
   1. run `trigeff` with `TRIGGER_SFS["2024"] = {"ele": None, "mu": "pog"}` (MUO reference for the
      muon trigger); electron cell -> `ele_trig_sf_2024.json.gz`
   2. set `"ele"` to that file, rerun `trigeff`; muon cell -> `muon_trig_sf_2024.json.gz`
   3. set `"mu"` to that file (`inference_common.TRIGGER_SFS`) for the analysis
   Until then `ele_trig_sf` = 1 and `mu_trig_sf` = MUO reference SF.
2. **b-tag SFs**: no shape SFs for UParT yet, so fixed-WP (method 1a) with our efficiency maps;
   `btag_sf` = 1 until then.
3. **Xbb working points**: no central GloParT WPs, so they come from the mistag rate of light-flavour
   fatjets (no B/C hadrons) in tt+LF events: 10% / 5% / 1% for LP / MP / HP.
   `run_pipeline.py inference --groups xbbwp`, then `python3 derive_xbb_wps.py --plots plots/xbbwp`
   prints the `analysis_years.XBB_WP_DICT` entry (and the ttH/ttZ(bb) efficiency, pT dependence).
   Then the multi-WP Xbb SFs as in Run 2.
4. **SPANet**: retrain on these definitions (`spanet_models.py`; inference runs without SPANet
   while `ACTIVE_MODELS` is empty), then pickling / datacards.

## SPANet training inputs

Everything that defines the training inputs is in `spanet_training.py`: samples and subsamples,
classes, the training category (`SELECTION`), AK4 cleaning, gen-matching dR, event weights, the
exported collections / fields (`EXPORT`) and the SPANet inputs (`SPANET_INPUTS`). After editing it:

```bash
python3 run_pipeline.py inference --backend pbs --groups spanet     # ntuples: outputs/nominal/spanet
python3 make_spanet_configs.py                                       # spanet/features_{parquet,h5}.yaml
cd ../spanet_inputs
python3 coffea_to_parquet_script.py -c ../Run3/spanet/features_parquet.yaml \
    -i ../Run3/outputs/nominal/spanet/output_all.coffea -o <dir>/run3.parquet --cat spanet
python3 parquet_to_h5.py -c ../Run3/spanet/features_h5.yaml -i <dir>/run3_*.parquet -o <dir>/run3.h5
```

- Workflow `spanet_training_workflow.py`: the analysis objects, then the training category and the gen
  matching (`parton_matching.py`). Partons come from the GenPart ancestry, not the LHE ordering, with
  the Run 2 provenance codes. AK4 jets carry `matched` / `prov`, and the fatjets are matched to the
  hard-process H/Z.
- `make_spanet_configs.py` warns if `SPANET_INPUTS` uses inputs the inference does not feed
  (`spanet_inference.onnx_inference`), and prints the `spanet_models.py` entry for the trained model.

## Running

```bash
python3 run_pipeline.py inference --backend pbs --groups trigeff
python3 run_pipeline.py inference --backend pbs --groups signal bkg minor data
python3 run_pipeline.py status    --backend pbs --groups trigeff signal bkg minor data
PC_SAMPLES=ttHTobb python3 run_pipeline.py inference --groups signal -- --executor futures --test
```
The PBS tools are shared with Run 2 (`../common`, see `../common/README.md`); the run options of
this analysis are in `run_options_pbs.yaml`.
