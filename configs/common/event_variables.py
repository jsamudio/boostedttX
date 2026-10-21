"""
Event-level SPANet inputs, one definition for Run 2 and Run 3, used by the training export
(spanet_training_workflow.py) and the inference (spanet_inference.onnx_inference):

  sphericity, aplanarity   from the eigenvalues l1 >= l2 >= l3 of the momentum tensor
                           S_ab = sum p_a p_b / sum |p|^2:  S = 3/2 (l2 + l3),  A = 3/2 l3
  centrality               sum pT / sum E
      over the leading lepton, the Z/H candidate (leading fat jet) and the AK4 jets outside it
      (dR > 0.8): the AK4 jets inside the fat jet are its constituents and would be counted twice
  n_jets                   number of AK4 jets
  n_bjets_<WP>             number of AK4 jets b-tagged at WP (JetGood.btag_<WP>, 0/1)
      all AK4 jets, or with count_clean only those outside the candidate

EVENT_FEATURES lists the names in the order the inputs are built; a model uses any subset.
"""
import awkward as ak
import numpy as np

CLEAN_DR = 0.8
SHAPE_FEATURES = ["sphericity", "aplanarity", "centrality"]


def event_features(btag_wps):
    """All available event-level feature names for the given b-tag working points."""
    return SHAPE_FEATURES + ["n_jets"] + [f"n_bjets_{wp}" for wp in btag_wps]


def outside_candidate(jets, fatjet):
    """AK4 jets with dR > 0.8 from the candidate (all jets if there is none)."""
    dr = jets.delta_r(fatjet)
    return jets[ak.fill_none(dr > CLEAN_DR, True)]


def _cartesian(objects):
    px = objects.pt * np.cos(objects.phi)
    py = objects.pt * np.sin(objects.phi)
    pz = objects.pt * np.sinh(objects.eta)
    energy = np.sqrt(px ** 2 + py ** 2 + pz ** 2 + objects.mass ** 2)
    return px, py, pz, energy


def event_shape(jets, fatjet, lepton):
    """(sphericity, aplanarity, centrality) as numpy arrays.
    jets: AK4 jets (jagged); fatjet, lepton: one object per event (may be None)."""
    def as_list(obj):
        return ak.singletons(obj)                      # None -> [], object -> [object]
    fields = ["pt", "eta", "phi", "mass"]
    parts = [ak.zip({f: x[f] for f in fields}) for x in
             (outside_candidate(jets, fatjet), as_list(fatjet), as_list(lepton))]
    objects = ak.concatenate(parts, axis=1)
    px, py, pz, energy = _cartesian(objects)
    p = [px, py, pz]
    norm = ak.to_numpy(ak.sum(px ** 2 + py ** 2 + pz ** 2, axis=1))
    tensor = np.zeros((len(norm), 3, 3))
    for a in range(3):
        for b in range(a, 3):
            tensor[:, a, b] = tensor[:, b, a] = ak.to_numpy(ak.sum(p[a] * p[b], axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        tensor /= norm[:, None, None]
        eig = np.sort(np.linalg.eigvalsh(np.nan_to_num(tensor)), axis=1)[:, ::-1]
        centrality = ak.to_numpy(ak.sum(objects.pt, axis=1)) / ak.to_numpy(ak.sum(energy, axis=1))
    sphericity = 1.5 * (eig[:, 1] + eig[:, 2])
    aplanarity = 1.5 * eig[:, 2]
    return sphericity, aplanarity, np.nan_to_num(centrality)


def event_variables(jets, fatjet, lepton, btag_wps, count_clean=False):
    """{feature name: numpy array} of event_features(btag_wps)."""
    sphericity, aplanarity, centrality = event_shape(jets, fatjet, lepton)
    counted = outside_candidate(jets, fatjet) if count_clean else jets
    out = {"sphericity": sphericity, "aplanarity": aplanarity, "centrality": centrality,
           "n_jets": ak.to_numpy(ak.num(counted, axis=1)).astype(np.float64)}
    for wp in btag_wps:
        field = f"btag_{wp}"
        if field not in counted.fields:
            raise ValueError(f"JetGood has no '{field}' (b-tag working points: {btag_wps})")
        out[f"n_bjets_{wp}"] = ak.to_numpy(ak.sum(counted[field], axis=1)).astype(np.float64)
    return out
