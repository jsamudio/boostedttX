#!/usr/bin/env python3
"""
SM weights of the 2024 tH samples from MiniAOD (run in a CMSSW environment: FWLite).

The THQ / THW ctcvcp samples are generated at kappa_t = -1 (the model's kHtt = +1) with the same 69-point
reweight card as the UL samples; the SM point (kHtt = -1, kSM = 1) is weight id rwgt_12 (index 11 of
LHEReweightingWeight in Run 2). NanoAODv15 dropped these weights (LHEReweightingWeight is empty), but
the MiniAOD LHEEventProduct has them (checked 2026-11-27 on one THQ file: 69 rwgt weights, mean of
rwgt_12 / nominal 0.0800, the minimum of the kSM = 1 scan).

Per MiniAOD file: (run, luminosityBlock, event) and w_SM = wgt(rwgt_12) / originalXWGTUP, the same
ratio as NanoAOD's LHEReweightingWeight. Written to th_sm_weights/<dataset>/<file>.npz, so a rerun skips
the files already done. --merge then writes th_sm_weights/<dataset>.npz (sorted keys, see
th_reweighting.py), which outputToPkl reads (the inference does not need it).

    source /cvmfs/cms.cern.ch/cmsset_default.sh && cmssw-el9 -B /cms/data,/home/$USER
    cmsenv (any CMSSW_15_0_X); voms-proxy-init -voms cms
    python3 extract_th_sm_weights.py --jobs 16                 # every file of THQ_2024 and THW_2024
    python3 extract_th_sm_weights.py --merge                   # then build the lookup tables
    python3 extract_th_sm_weights.py --max-files 2             # quick test first

Only the LHEEventProduct branch is read (remotely, --redirector), so a file costs much less than its size.
"""
import argparse
import glob
import multiprocessing as mp
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "th_sm_weights")
SM_WEIGHT_ID = "rwgt_12"         # kHtt = -1, kSM = 1 (SM); index 11 of the 69 points
N_POINTS = 69
MINIAOD = {
    "THQ_2024": "/THQ-4FS-ctcvcp_Par-M-125_TuneCP5_13p6TeV_madgraph-pythia8/"
                "RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM",
    "THW_2024": "/THW-5FS-ctcvcp_Par-M-125_TuneCP5_13p6TeV_madgraph-pythia8/"
                "RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM",
}


def das_files(das_name):
    out = subprocess.run(["dasgoclient", "-query", f"file dataset={das_name}"],
                         capture_output=True, text=True, check=True).stdout.split()
    if not out:
        raise SystemExit(f"No files in DAS for {das_name}")
    return sorted(out)


def file_out(dataset, lfn):
    return os.path.join(OUT_DIR, dataset, os.path.splitext(os.path.basename(lfn))[0] + ".npz")


def extract(task):
    dataset, lfn, redirector = task
    out = file_out(dataset, lfn)
    if os.path.exists(out):
        return dataset, lfn, "done before"
    from DataFormats.FWLite import Events, Handle          # imported in the worker (spawn)
    handle = Handle("LHEEventProduct")
    runs, lumis, evts, ws, signs = [], [], [], [], []
    pos = None
    for e in Events(redirector + lfn):
        e.getByLabel("externalLHEProducer", handle)
        p = handle.product()
        weights = p.weights()
        if pos is None:
            ids = [weights[i].id for i in range(weights.size())]
            rwgt = [i for i, x in enumerate(ids) if x.startswith("rwgt_")]
            if len(rwgt) != N_POINTS or SM_WEIGHT_ID not in ids:
                raise RuntimeError(f"{lfn}: {len(rwgt)} rwgt weights (expected {N_POINTS}), "
                                   f"{SM_WEIGHT_ID} {'present' if SM_WEIGHT_ID in ids else 'missing'}")
            pos = ids.index(SM_WEIGHT_ID)
        w = weights[pos]
        if w.id != SM_WEIGHT_ID:                     # same layout in every event of a file, checked
            raise RuntimeError(f"{lfn}: weight {pos} is {w.id}, not {SM_WEIGHT_ID}")
        aux = e.eventAuxiliary()
        runs.append(aux.run())
        lumis.append(aux.luminosityBlock())
        evts.append(aux.event())
        ws.append(w.wgt / p.originalXWGTUP())
        signs.append(1 if p.originalXWGTUP() >= 0 else -1)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".tmp.npz"
    np.savez(tmp, run=np.asarray(runs, dtype=np.uint32), lumi=np.asarray(lumis, dtype=np.uint32),
             event=np.asarray(evts, dtype=np.uint64), w=np.asarray(ws, dtype=np.float32),
             sign=np.asarray(signs, dtype=np.int8))
    os.replace(tmp, out)
    return dataset, lfn, f"{len(ws)} events, mean w_SM {np.mean(ws) if ws else float('nan'):.4f}"


def merge(dataset):
    """th_sm_weights/<dataset>.npz: key = (luminosityBlock << 32) | event, sorted, with w and the runs."""
    from th_reweighting import event_keys
    parts = sorted(glob.glob(os.path.join(OUT_DIR, dataset, "*.npz")))
    if not parts:
        raise SystemExit(f"No extracted files for {dataset} in {OUT_DIR}/{dataset}")
    run = np.concatenate([np.load(p)["run"] for p in parts])
    lumi = np.concatenate([np.load(p)["lumi"] for p in parts])
    event = np.concatenate([np.load(p)["event"] for p in parts])
    w = np.concatenate([np.load(p)["w"] for p in parts])
    sign = np.concatenate([np.load(p)["sign"] for p in parts])
    keys = event_keys(lumi, event)
    order = np.argsort(keys, kind="stable")
    keys, w, sign = keys[order], w[order], sign[order]
    dup = np.count_nonzero(keys[1:] == keys[:-1])
    if dup:
        raise SystemExit(f"{dataset}: {dup} duplicate (lumi, event) keys")
    path = os.path.join(OUT_DIR, f"{dataset}.npz")
    np.savez(path, key=keys, w=w, sign=sign, runs=np.unique(run))
    print(f"{dataset}: {len(parts)} files, {len(keys)} events, runs {np.unique(run).tolist()}, "
          f"mean w_SM {w.mean():.4f} -> {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", nargs="+", default=list(MINIAOD), choices=list(MINIAOD))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--redirector", default="root://cms-xrd-global.cern.ch/")
    ap.add_argument("--max-files", type=int, default=None, help="first N files per dataset (test)")
    ap.add_argument("--merge", action="store_true", help="build the lookup tables from the extracted files")
    args = ap.parse_args()

    if args.merge:
        for ds in args.datasets:
            merge(ds)
        return
    tasks = []
    for ds in args.datasets:
        files = das_files(MINIAOD[ds])[:args.max_files]
        todo = [f for f in files if not os.path.exists(file_out(ds, f))]
        print(f"{ds}: {len(files)} files, {len(todo)} to do")
        tasks += [(ds, f, args.redirector) for f in todo]
    failed = []
    with mp.get_context("spawn").Pool(args.jobs) as pool:
        for res in pool.imap_unordered(_safe_extract, tasks):
            print(" ", *res)
            if res[2].startswith("FAILED"):
                failed.append(res[1])
    if failed:
        print(f"{len(failed)} files failed; rerun to retry them (done files are skipped)")


def _safe_extract(task):
    try:
        return extract(task)
    except Exception as err:                         # noqa: BLE001  (report, keep the other files going)
        return task[0], task[1], f"FAILED: {err}"


if __name__ == "__main__":
    main()
