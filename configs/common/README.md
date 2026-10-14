# Shared pocket-coffea tools

Analysis-independent pieces used by the analyses in `configs/` (requires the latest
pocket-coffea, 0.9.x with the calibrator system).

| File | What |
|---|---|
| `executors_pbs.py` | pocket-coffea executor for PBS Pro: one batch job per chunk of files, each running the futures executor on all cores of its node inside apptainer |
| `pbs_tools.py` | job array + merge job submission, status, resubmission, VOMS proxy refresh (stdlib only, runs on the host) |
| `run_options_pbs.yaml` | template run options (queue, cores, memory, image, binds, splitting) |
| `met_calibration.py` | type-1 PF MET for Run 2 UL, consistent in data and MC (replaces pocket-coffea's `METCalibrator`) |

## Running an analysis on PBS

Prerequisites: a recent pocket-coffea config; the analysis' own modules live in its config
directory; input files readable from the compute nodes (under the apptainer binds).

1. Copy `run_options_pbs.yaml` into your config directory and adapt it (mainly
   `max-events-per-job`, `chunksize`, `mem-per-worker`).
2. VOMS proxy (inside the apptainer image if `voms-proxy-init` is not on the host):
   ```bash
   voms-proxy-init -voms cms -rfc --valid 168:0
   ```
3. Prepare the jobs, from your config directory inside the apptainer image:
   ```bash
   pocket-coffea run --cfg <config>.py -o <outputdir> --executor pbs \
       --executor-custom-setup ../common/executors_pbs.py -ro run_options_pbs.yaml
   ```
4. Submit and follow them from the host (where `qsub` is):
   ```bash
   python3 ../common/pbs_tools.py submit   <outputdir>
   python3 ../common/pbs_tools.py status   <outputdir>
   python3 ../common/pbs_tools.py resubmit <outputdir>     # failed jobs + a new merge
   ```
   The merge job writes `<outputdir>/output_all.coffea`; job logs are in `<outputdir>/job/logs/`.

What a resubmission picks up: jobs re-import the modules of the analysis directory and of
`configs/common` from disk, but the configuration (and modules registered with
`cloudpickle.register_pickle_by_value`) is frozen when the jobs are prepared. For changes there,
delete `<outputdir>` and prepare again.

## Type-1 MET (`met_calibration.py`)

pocket-coffea's default `METCalibrator` rebuilds MET from `RawMET`: in data (`apply_data: false`)
that is the raw MET, and in MC it uses no L1 offset, which is wrong for Run 2 CHS jets. To keep
NanoAOD's type-1 MET and only propagate the changes of the calibrators (AK4 JEC/JER/JES, lepton
scale), pass the sequence from this module to the Configurator:
```python
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
from met_calibration import CALIBRATORS
cfg = Configurator(..., calibrators=CALIBRATORS, ...)
```

## Notes for new analyses

* Data from several primary datasets: give each data sample a trigger requirement that removes
  the events already contained in the other datasets (see `hlt_selection` in
  `Run2/custom_cut_functions.py` for the single-lepton case), otherwise overlapping events are
  counted more than once.
* Per-sample (`bysample`) columns/weights/variations work with the PBS executor: each job's
  configuration is restricted to the samples it processes.
