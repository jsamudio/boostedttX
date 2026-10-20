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
   PC_JEC_VARIATION=JER pocket-coffea run --cfg inference_conf_minor.py -o outputs/JER/minor
   PC_SAMPLES=ttHTobb pocket-coffea run --cfg inference_conf_signal.py -o test --test   # quick test
   ```
3. **b-tag SF normalization**: `python btag_norm.py -i outputs/nominal/*/output_all.coffea -o btag_norm.json`
4. **Pickles**: `python outputToPkl.py -i outputs/nominal/*/output_all.coffea --btag-norm btag_norm.json`
   (jet-variation outputs: add `--variations all --skip-nominal`)
5. **Datacards**: `python makeDatacard.py [--eft-ttbb]`

`run_pipeline.py` chains steps 2-5 with the output layout above:
```bash
python run_pipeline.py inference --groups signal bkg minor --eras 2018 --jec-variations JER
python run_pipeline.py btagnorm pickle datacard --jec-variations JER
python run_pipeline.py all --dry-run            # show the commands only
python run_pipeline.py inference -- --executor futures --test   # args after -- go to pocket-coffea
```
(see below for running the inference as PBS jobs)

## Running on the PBS Pro cluster

Each group / jet variation becomes a PBS job array: pocket-coffea splits the files into jobs
(`max-events-per-job` in `run_options_pbs.yaml`), every array subjob runs the futures executor on
8 cores inside the pocketcoffea apptainer image, and a merge job (queued after the array)
writes `output_all.coffea`. Settings (queue, cores, memory, image, binds) are in `run_options_pbs.yaml`.

The jobs need a VOMS proxy. Create it first (inside the apptainer image if `voms-proxy-init` is
not available on the host; it lands in `/tmp/x509up_u<uid>`, visible to both):
```bash
voms-proxy-init -voms cms -rfc --valid 168:0
```
Every `submit` / `resubmit` copies the current proxy to `$HOME` (exported to the jobs as
`X509_USER_PROXY`) and refuses to submit if it has less than 1 h left, so renew it before
resubmitting.

From the host (where `qsub` is available; pocket-coffea is run inside apptainer automatically):
```bash
python run_pipeline.py inference --backend pbs --groups signal bkg minor --eras 2018 --jec-variations JES_ALL JER JMSJMR
python run_pipeline.py status    --backend pbs --jec-variations JES_ALL JER JMSJMR     # done/failed/running per group
python run_pipeline.py resubmit  --backend pbs --jec-variations JES_ALL JER JMSJMR     # failed jobs + new merge
python run_pipeline.py btagnorm pickle datacard --jec-variations JES_ALL JER JMSJMR    # once everything is merged
```
`resubmit --include-unfinished` also resubmits jobs that never finished (e.g. killed by the batch
system); it refuses while PBS still lists the arrays. Per directory the same is available as
`python ../common/pbs_tools.py submit|status|resubmit|merge <outputdir>`. Job logs: `<outputdir>/job/logs/job_<id>.log`.
The PBS executor and tools are shared with other analyses in `../common` (see `../common/README.md`);
`run_options_pbs.yaml` here holds this analysis' settings.

## pocket-coffea version

The configs target the latest pocket-coffea (0.9.x, calibrator-based). `params/*.yaml` only hold
overrides of the pocket-coffea defaults (correction files come from cvmfs via `${cvmfs:...}`):
muon SF names for the MediumID selection, the UL jet veto maps and the MET filters.
Shape variations are selected per job with `PC_JEC_VARIATION` (`inference_common.py`):
`JES_<source>` (each individual JES source of `analysis_years.JES_SOURCES`, AK4 + AK8 jets shifted
together, checked against the JME files), `JER` (AK4 and AK8 varied separately) and `JMSJMR`
(AK8 mass scale/resolution, declared by `fatjet_mass_variations.py`, applied in the processor). MET filters and the golden JSON are applied with the
standard `eventFlags` / `goldenJson` skim cuts; JEC/JER and type-1 MET come from the default
calibrator sequence.

## Notes on the inference workflow

* The b-tag requirement (`btag_mask`, >= 2 b jets) is a category, not a preselection: the
  b-tag SF normalization needs the events before it. The processor accumulates the
  normalization sums per jet-multiplicity bin in `output['btag_calib']` (no per-event
  columns), then drops the events failing `btag_mask` before the expensive part
  (ZH candidate variables, control variables, SPANet). `btag_norm.py` reads these sums;
  older outputs with a `btag_calib` column category still work.
* For JES shape variations (e.g. `AK4AK8_JES_AbsoluteStatUp`) `btag_sf` is the JES-varied BTV
  shape SF (`up_jes*`/`down_jes*`, c jets keep the central SF), and `btag_norm.py` derives the
  normalization separately for every shape variation (`shape_variations` in the JSON).
* Nominal MC weight (`outvars.nominal_weight_factors`) includes `l1prefire_sf` (UL NanoAOD
  `L1PreFiringWeight`, ECAL + muon prefiring) for all years; `l1prefire_<year>` is a datacard shape systematic.
* The 2018 HEM15/16 veto is in the preselection of every group (no effect on other eras).
* JECs: the recommended JECs are reapplied to data and MC jets (`params/jets_calibration.yaml`,
  `apply_jec_Data`), as JME requires.
* All MET used (selection, SPANet inputs, `MET_pt`/`MET_phi` columns) is the type-1 PF `MET`,
  recomputed from `RawMET` with the same corrected AK4 jets (`Jet`, `CorrT1METJet`) and the L1 of the
  same payload, for every jet variation, by `../common/met_calibration.py` (replaces pocket-coffea's
  `METCalibrator`, which uses no L1 for CHS jets and returns raw MET in data).
* Fatjet mass: every fatjet mass is the ParticleNet regressed mass, JMS/JMR-corrected in MC and as-is
  in data (`custom_scale_factors.set_fatjet_mass`): `FatJet.mass` and `FatJet.particleNet_mass` both
  carry it, so the FatJetGood mass window, `ZH_M`, the fatjet mass columns, the ZH candidate
  4-vector and the Xbb SF binning all use the same value.
* Triggers (`custom_cut_functions.hlt_selection`): MC uses the OR of the e/gamma and muon paths;
  data uses the paths of its primary dataset, with SingleMuon events removed from the e/gamma
  dataset, so no data event is counted twice.
* Correction files, efficiency maps and ONNX sessions are cached per dask worker.
* The JMR stochastic smearing is seeded per chunk, so it is reproducible and identical
  for all shape variations.
