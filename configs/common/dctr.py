"""
DCTR reweighting of the ttbar hdamp variation: per-event weights that turn the nominal ttbar sample
into the hdamp up / down prediction, from two classifiers (ONNX) trained on the gen top quarks.

Inputs per event: the gen top and antitop (GenPart, pdgId +-6, isFirstCopy), six features each
    log10(pT), rapidity, phi, mass / MAX_MASS, particle id (t 0.1, tbar 0.2), hdamp / mt (nominal)
stacked as (t, tbar) -> shape (n, 2, 6). The weight is p0 / p1 of the classifier output.
Events with pT(ttbar) >= MAX_TTBAR_PT, or without exactly one top and one antitop, get weight 1.

Optional: the processors call it only for the samples given in the workflow option 'dctr_hdamp'
(inference_common.py, PC_DCTR=1), before the skim, so that the sums over ALL generated events
can normalize the variations to the nominal cross section:

    output["dctr_sums"][sample][dataset][region][variation] = sum of sign(genWeight) * weight

  variation : nominal (weight 1), dctr_hdamp_up, dctr_hdamp_down
  region    : all, and ttB (genTtbarId % 100 >= 51) for samples with genTtbarId
"""
import awkward as ak
import numpy as np
import vector

vector.register_awkward()

MAX_MASS = 243.95          # top mass normalization of the training
HDAMP_NOMINAL = 1.379      # hdamp / mt of the nominal samples (237.9 GeV / 172.5 GeV)
MAX_TTBAR_PT = 1000.       # GeV, the classifiers are not used above
VARIATIONS = ("dctr_hdamp_up", "dctr_hdamp_down")

_LOCAL_CACHE = {}


def _session(path):
    """ONNX session, kept on the dask worker (or in the module) between chunks."""
    try:
        from dask.distributed import get_worker
        worker = get_worker()
        if not hasattr(worker, "dctr_cache"):
            worker.dctr_cache = {}
        cache = worker.dctr_cache
    except (ImportError, ValueError):
        cache = _LOCAL_CACHE
    if path not in cache:
        import onnxruntime as ort
        opts = ort.SessionOptions()
        opts.intra_op_num_threads = 1
        opts.inter_op_num_threads = 1
        opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        cache[path] = ort.InferenceSession(path, sess_options=opts, providers=['CPUExecutionProvider'])
    return cache[path]


def _features(top, pid):
    n = len(top)
    return np.stack([np.log10(ak.to_numpy(top.pt)), ak.to_numpy(top.rapidity), ak.to_numpy(top.phi),
                     ak.to_numpy(top.mass) / MAX_MASS, np.full(n, pid), np.full(n, HDAMP_NOMINAL)], axis=-1)


def _ratio(path, inputs):
    sess = _session(path)
    pred = sess.run([sess.get_outputs()[0].name], {sess.get_inputs()[0].name: inputs})[0]
    return pred[:, 0] / pred[:, 1]


def dctr_hdamp_weights(events, models):
    """(up, down) per-event weights; models = {"up": onnx path, "down": onnx path}."""
    up, down = np.ones(len(events)), np.ones(len(events))
    gen = events.GenPart
    first = gen.hasFlags(['isFirstCopy'])
    tops, antitops = gen[(gen.pdgId == 6) & first], gen[(gen.pdgId == -6) & first]
    valid = ak.to_numpy((ak.num(tops) == 1) & (ak.num(antitops) == 1))
    if not valid.any():
        return up, down

    def p4(p):
        return vector.zip({"pt": p.pt, "eta": p.eta, "phi": p.phi, "mass": p.mass})
    t, tbar = p4(tops[valid][:, 0]), p4(antitops[valid][:, 0])
    inputs = np.stack([_features(t, 0.1), _features(tbar, 0.2)], axis=1).astype(np.float32)
    in_range = ak.to_numpy((t + tbar).pt) < MAX_TTBAR_PT
    up[valid] = np.where(in_range, _ratio(models["up"], inputs), 1.0)
    down[valid] = np.where(in_range, _ratio(models["down"], inputs), 1.0)
    return up, down


def fill_dctr_sums(events, weights, output, sample, dataset):
    """Add this chunk's sums to output['dctr_sums'][sample][dataset][region][variation]."""
    sign = np.sign(ak.to_numpy(events.genWeight))
    regions = {"all": np.ones(len(events), dtype=bool)}
    if "genTtbarId" in events.fields:
        regions["ttB"] = ak.to_numpy(events.genTtbarId % 100 >= 51)
    out = output.setdefault("dctr_sums", {}).setdefault(sample, {}).setdefault(dataset, {})
    for region, mask in regions.items():
        sums = out.setdefault(region, {})
        for variation, w in {"nominal": np.ones(len(events)), **weights}.items():
            sums[variation] = sums.get(variation, 0.0) + float(np.sum(sign[mask] * w[mask]))


def apply_dctr_hdamp(events, output, sample, dataset, options):
    """Processor hook (before the skim). options = workflow option 'dctr_hdamp' or None.
    Sets events['dctr_hdamp_up'/'dctr_hdamp_down'] and the sums for the samples in options['samples']."""
    if not options or sample not in options["samples"] or len(events) == 0:
        return
    weights = dict(zip(VARIATIONS, dctr_hdamp_weights(events, options["models"])))
    for name, w in weights.items():
        events[name] = w
    fill_dctr_sums(events, weights, output, sample, dataset)


def set_default_dctr_columns(events):
    """Weight 1 for the samples without DCTR weights, so that the exported columns exist for all MC."""
    for name in VARIATIONS:
        if name not in events.fields:
            events[name] = np.ones(len(events))
