# Run 3 semileptonic boosted ttZ/ttH analysis (2024)

Same structure as `configs/Run2`, with the Run 3 object definitions. The previous 2024 setup
(`configs/SemiLeptonic_2024`) is kept unchanged for comparison until this one is validated.

## Object definitions

| | 2024 |
|---|---|
| Data / MC | NanoAODv15, Summer24 MC, runs C-I (109.82 fb^-1, golden JSON) |
| MET | PUPPI MET; pocket-coffea's MET calibrator recomputes the type-1 PUPPI MET from `RawPuppiMET` with the recommended JECs (data and MC) |
| JECs | pocket-coffea defaults: `Summer24Prompt24_V5`, JER `Summer24Prompt24_JRV2`, reapplied to data and MC (AK4 PUPPI, CorrT1METJet, AK8) |
| AK4 b tagging | UParTAK4, medium WP from pocket-coffea's btagging parameters (`object_preselection.btag.wp`); fixed-WP SFs over L/M/T (multi-WP method 1a) |
| AK8 tagger | GloParT XbbVsQCD = Xbb / (Xbb + QCD); the ZH candidate is the leading fatjet in it |
| Fatjet mass | GloParT 2-prong regressed mass (`globalParT3_massCorrX2p` x raw mass), in data and MC, for the mass window, `ZH_M`, the columns and 4-vectors (`custom_scale_factors.set_fatjet_mass`); no nominal JMS/JMR (JMAR, DNN-regressed masses), JMS 2% / JMR 3% as the `JMSJMR` variations (`--jec-variations JMSJMR`, `../common/fatjet_mass_variations.py`) |
| Triggers | Ele30_WPTight_Gsf, Ele115_CaloIdVT_GsfTrkIdT / IsoMu24, Mu50, HighPtTkMu100, CascadeMu100; per primary dataset in data (Muon PD: muon paths; EGamma PD: e/gamma paths and no muon path) |
| Leptons | electrons MVA noIso WP90, muons medium ID, mini-isolation (`params/object_preselection.yaml`) |

Process labels come from `samples.PROCESS_LABELS` (ttZ = TTZToQQ incl. Z->bb, TTLL, TTNuNu;
ttH = ttHTobb, ttHToNonbb); ttbar samples are split by genTtbarId.

## Not derived yet (bottom-up)

1. **Trigger SFs** (next): ttbar e-mu without trigger requirement + `TrigEff.ipynb` (orthogonal
   method), two groups as in Run 2:
   1. `trigeff_ele` (MC + Muon data): the reference muon is the MUO trigger-SF denominator (medium ID +
      PF isolation medium, `pfIsoId >= 3`, no sip3d) and gets the MUO muon trigger SF, whatever
      `TRIGGER_SFS` holds; electron cell (`COFFEA_DIR = outputs/nominal/trigeff_ele/`) ->
      `ele_trig_sf_2024.json.gz`
   2. set `TRIGGER_SFS["2024"]["ele"]` to that file; `trigeff_mu` (MC + EGamma data, analysis muon:
      mini-isolation, sip3d); muon cell (`outputs/nominal/trigeff_mu/`) -> `muon_trig_sf_2024.json.gz`
   3. set `"mu"` to that file (`inference_common.TRIGGER_SFS`) for the analysis
   Until then `ele_trig_sf` = 1 and `mu_trig_sf` = MUO reference SF.

       python3 run_pipeline.py inference --backend pbs --groups trigeff_ele
       python3 run_pipeline.py inference --backend pbs --groups trigeff_mu    # after step 2
2. **b-tag SFs**: fixed-WP SFs (BTV `UParTAK4_comb` for b/c, `UParTAK4_light`), multi-WP method 1a
   over the WPs whose information the analysis uses (`analysis_years.BTAG_SF_WPS` = L/M/T: the
   `n_bjets_<WP>` counts and the M requirement), on all `JetGood`, with our MC efficiency maps
   (`custom_scale_factors.sf_btag_multiwp`; a factor with a negative numerator is set to 1, as BTV
   recommends). Maps per process group (`samples.BTAG_EFF_GROUPS`):
   1. `run_pipeline.py inference --groups btageff` (preselection, before the b-tag requirement)
   2. `python3 derive_btag_eff.py --plots plots/btageff` -> `btag_eff_2024.json.gz`
   3. set `inference_common.BTAG_SFS["2024"]["eff"]` to that file
   Columns: `btag_sf` and `btag_sf_{bc,light}_{correlated,uncorrelated}_{up,down}`; `btag_sf` = 1
   until the maps exist. SPANet also uses the continuous `btagB`, which these SFs do not correct.
3. **Xbb working points**: no central GloParT WPs, so they come from the mistag rate of light-flavour
   fatjets (no B/C hadrons) in tt+LF events: 10% / 5% / 1% for LP / MP / HP.
   `run_pipeline.py inference --groups xbbwp`, then `python3 derive_xbb_wps.py --plots plots/xbbwp`
   prints the `analysis_years.XBB_WP_DICT` entry (and the ttH/ttZ(bb) efficiency, pT dependence).
   On the 2024 outputs the light-mistag targets give very loose cuts (LP 0.0026, 99.7% signal
   efficiency), so the WPs are to be chosen by sensitivity instead: `scan_xbb_wp.py` reads the main
   inference output (signal, bkg, minor) and prints the expected sigma(mu) of ttZ / ttH / tt_B (Fisher,
   pT x mass bins) per cut (`cut` mode), and ranks sets of WP categories above a preselection against
   the continuous score (`multiwp` mode, for the SPANet multi-WP Xbb input).
   Then the multi-WP Xbb SFs as in Run 2.
4. **SPANet**: retrained on these definitions (`spanet_models.py`: `run3_refactor`, ONNX in
   `spanet/models/`; the inference runs without SPANet if `ACTIVE_MODELS` is empty), then pickling /
   datacards.
5. **Pickling and datacards**: ported from Run 2 (see below); the 13.6 TeV inputs listed there are
   still placeholders.

## Pickling and datacards (ported from Run 2)

Same chain as Run 2 without `btag_norm.py` (the fixed-WP SFs need no normalization factor):

    python3 run_pipeline.py pickle norms datacard --groups signal bkg minor [--jec-variations JMSJMR ...]

| Script | Run 3 differences from Run 2 |
|---|---|
| `outputToPkl.py` | 13.6 TeV `XSEC_DICT` by dataset name (variation samples: the nominal dataset, `VARIATION_NOMINAL`); no `--btag-norm`; no tH SM reweighting (2024 THQ / THW taken as SM); `VV` process |
| `make_process_norms.py`, `theory_norms.py` | ttZ = TTZToQQ (with its Z -> bb events) + TTLL + TTNuNu; ttX sample `TTX`. The processor fills `theory_sums` before the skim (`process_extra_before_skim`) |
| `makeDatacard.py` | VV background; b-tag nuisances `btag_{bc,light}_correlated` and `btag_{bc,light}_uncorrelated_<year>`; muon ID / iso total per year (no stat / syst split, no reco SF); no Xbb SF, L1 prefiring or PU jet ID; jet nuisances found from the varied pickles (`jet_systematics`: JES correlated, `jes_<source>_ak4` / `_ak8` separately since pocket-coffea varies AK4 and AK8 jets as separate variations; JER, JMS, JMR per year) |
| `datacard_shapes.py`, `pdf_uncertainty.py`, `TH1.py`, `plot_nn_edges.py` | copies of the Run 2 files |

The processor writes for this: `pdfweight_0..102` (signal and tt+jets), `tt_2b` / `tt_C`
(tt+2b / tt+C cross-section nuisances) and `theory_sums`.

Placeholders (printed by `makeDatacard.py` at every run, `PLACEHOLDERS`): `THEORY_LNN` (signal and
ttbar cross-section uncertainties, 13 TeV values), `TTBAR_MODELLING_LNN` (UE, hdamp, hdamp_ttbb),
`BKG_XSEC_UNC` (VV, single top), `VJETS_*_UNC`, `analysis_years.LUMI_UNCERTAINTIES` (empty: no lumi
nuisance yet). `outputToPkl.XSEC_DICT` holds the values of the dataset definitions, to be checked
against their references. `pdf_uncertainty.py` knows the UL PDF sets (103 / 101 members); a Summer24
sample with another set stops the datacard step there. `analysis_years.XBB_WP_DICT` (the LP is the
datacard cut) is provisional: LP 0.25, M / T as the SPANet model.

## SPANet training inputs

Everything that defines the training inputs is in `spanet_training.py`: samples and subsamples,
classes (index = SPANet output position, 0..N-1, as in `spanet_models.SPANET_CLASSES`), the training
category (`SELECTION`), AK4 cleaning, the event-level inputs, gen-matching dR, event weights, the
exported collections / fields (`EXPORT`) and the SPANet inputs (`SPANET_INPUTS`). After editing it:

```bash
python3 run_pipeline.py inference --backend pbs --groups spanet     # ntuples: outputs/nominal/spanet
python3 ../common/make_spanet_configs.py                             # spanet/*.yaml
cd ../spanet_inputs
python3 coffea_to_parquet_script.py -c ../Run3/spanet/features_parquet.yaml \
    -i ../Run3/outputs/nominal/spanet/output_all.coffea -o <dir>/run3.parquet --cat spanet
python3 parquet_to_h5.py -c ../Run3/spanet/features_h5.yaml -i <dir>/run3_*.parquet -o <dir>/run3.h5
```

- Workflow `spanet_training_workflow.py`: the analysis objects, the selection variables (all AK4
  jets, as the analysis), the event-level inputs, then the training category and the gen matching
  (`../common/parton_matching.py`, shared with Run 2). Partons come from the GenPart ancestry, not the
  LHE ordering. AK4 jets carry `matched` / `prov`, and the fatjets are matched to the hard-process H/Z.
- Event-level inputs (`../common/event_variables.py`, SPANet input `Event`): sphericity, aplanarity
  and centrality over the lepton, the candidate fat jet and the AK4 jets outside it; the number of AK4
  jets and of b-tagged jets per working point (`COUNT_CLEAN_JETS`: only jets outside the candidate).
- `make_spanet_configs.py` also writes `spanet/event_info_inputs.yaml` (the INPUTS section of the
  SPANet event-info file), warns if `SPANET_INPUTS` uses inputs the inference does not feed, and prints
  the `spanet_models.py` entry for the trained model.
- Aligned with the Run 2 model `refactor_noOverlap_MultiWP` (2026-11-24): signal = gen-matched
  `ttHTobb` / `TTZToQQ` (`matchedGen_ZHbb_bb`, the Run 2 flag; it also selects Z -> bb in TTZToQQ), AK4
  inputs only outside the candidate (`CLEAN_JETS = True`), fat-jet inputs `xbb_wp_L/M/T`, training cut
  candidate >= L. Intended differences: UParT, PUPPI MET, GloParT, no TTZ->bb sample.
- Training: `make_spanet_configs.py` also writes the complete event-info file `spanet/event_info.yaml`
  (`EVENT_TOPOLOGY`); options `../spanet_inputs/options_files/options_run3.json` (Run 2 hyperparameters,
  no file paths: SPANet's options file would override the command line). Validation = last 5% of the
  training h5 (`train_validation_split`), so the converter's test h5 stays independent for
  `../spanet_inputs/overtraining_check.py` (train vs test and validation vs test, per class).
- Xbb working points are not in the export: the parquet keeps the continuous `FatJetSorted.bbtag` and no
  Xbb cut (`SELECTION["min_xbb"] = None`). The training's WPs (`XBB_WPS`, inputs `xbb_wp_L/M/T` =
  bbtag >= cut) and candidate cut (`TRAIN_MIN_XBB`) are applied by `parquet_to_h5.py` (section `xbb` of
  `features_h5.yaml`). Another WP set = edit them, rerun `make_spanet_configs.py` and `parquet_to_h5.py`
  on the same parquet. The model entry records its `xbb_wps`; the inference builds the flags per model.

## Running

```bash
python3 run_pipeline.py inference --backend pbs --groups trigeff_ele btageff
python3 run_pipeline.py inference --backend pbs --groups signal bkg minor data
python3 run_pipeline.py status    --backend pbs --groups trigeff_ele signal bkg minor data
PC_SAMPLES=ttHTobb python3 run_pipeline.py inference --groups signal -- --executor futures --test
```
The PBS tools are shared with Run 2 (`../common`, see `../common/README.md`); the run options of
this analysis are in `run_options_pbs.yaml`.
