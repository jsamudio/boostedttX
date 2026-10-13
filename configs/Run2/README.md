# Run 2 semileptonic boosted ttZ/ttH analysis

## Where things are configured

| What | File |
|---|---|
| Eras, luminosity, Xbb working points, pickle directory, years to make datacards for | `analysis_years.py` |
| Samples, sample groups, subsample splits, which subsamples are pickled | `samples.py` |
| SPANet ONNX models and which ones are active | `spanet_models.py` |
| Inference config (selection, columns, calibrators, workflow options) | `inference_common.py` |
| Inference workflow | `spanet_inference.py` |
| Output columns, weight factors | `outvars.py` |

`inference_conf_{signal,bkg,minor,EFT,data}.py` are two-line wrappers around
`inference_common.build_config(<group>)`.

## Workflow

1. **SPANet training inputs**: `pocket-coffea run --cfg testconf_spanet.py -o spanet_inputs`
   (workflow `zhbbworkflow_spanet.py`), then train and export the ONNX model and add it to
   `spanet_models.py` (set `clean_jets` to match the training inputs).
2. **Inference**, one run per group:
   ```bash
   pocket-coffea run --cfg inference_conf_signal.py -o outputs/nominal/signal
   PC_ERAS=2017,2018 pocket-coffea run --cfg inference_conf_bkg.py -o outputs/nominal/bkg
   PC_JEC_VARIATION=JES_Total pocket-coffea run --cfg inference_conf_minor.py -o outputs/JES_Total/minor
   PC_SAMPLES=ttHTobb pocket-coffea run --cfg inference_conf_signal.py -o test --test   # quick test
   ```
3. **b-tag SF normalization**: `python btag_norm.py -i outputs/nominal/*/output_all.coffea -o btag_norm.json`
4. **Pickles**: `python outputToPkl.py -i outputs/nominal/*/output_all.coffea --btag-norm btag_norm.json`
   (jet-variation outputs: add `--variations all --skip-nominal`)
5. **Datacards**: `python makeDatacard.py [--eft-ttbb]`

`run_pipeline.py` chains steps 2-5 with the output layout above:
```bash
python run_pipeline.py inference --groups signal bkg minor --eras 2018 --jec-variations JES_Total
python run_pipeline.py btagnorm pickle datacard --jec-variations JES_Total
python run_pipeline.py all --dry-run            # show the commands only
python run_pipeline.py inference -- --executor futures --test   # args after -- go to pocket-coffea
```

## pocket-coffea version

The configs target the latest pocket-coffea (0.9.x, calibrator-based). `params/*.yaml` only hold
overrides of the pocket-coffea defaults (correction files come from cvmfs via `${cvmfs:...}`):
muon SF names for the MediumID selection, the UL jet veto maps, the MET filters, and the jet
variations (selected per job with `PC_JEC_VARIATION`; names are checked against
`default_jets_calibration.variations`). MET filters and the golden JSON are applied with the
standard `eventFlags` / `goldenJson` skim cuts; JEC/JER and type-1 MET come from the default
calibrator sequence.

## Notes on the inference workflow

* The b-tag requirement (`btag_mask`, >= 2 b jets) is a category, not a preselection: the
  b-tag SF normalization needs the events before it. The processor accumulates the
  normalization sums per jet-multiplicity bin in `output['btag_calib']` (no per-event
  columns), then drops the events failing `btag_mask` before the expensive part
  (ZH candidate variables, control variables, SPANet). `btag_norm.py` reads these sums;
  older outputs with a `btag_calib` column category still work.
* For AK4 JES shape variations (e.g. `AK4PFchs_JES_TotalUp`) `btag_sf` is the JES-varied BTV
  shape SF (`up_jes*`/`down_jes*`, c jets keep the central SF), and `btag_norm.py` derives the
  normalization separately for every shape variation (`shape_variations` in the JSON).
* Nominal MC weight (`outvars.nominal_weight_factors`) includes `l1prefire_sf` (UL NanoAOD
  `L1PreFiringWeight`, ECAL + muon prefiring) for all years; `l1prefire_<year>` is a datacard shape systematic.
* The 2018 HEM15/16 veto is in the preselection of every group (no effect on other eras).
* All MET used (selection, SPANet inputs, `MET_pt`/`MET_phi` columns) is the type-1 corrected PF `MET`.
* Correction files, efficiency maps and ONNX sessions are cached per dask worker.
* The JMR stochastic smearing is seeded per chunk, so it is reproducible and identical
  for all shape variations.
