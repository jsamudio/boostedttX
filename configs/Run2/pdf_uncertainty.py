"""
PDF and alpha_s uncertainties from the LHEPdfWeight members. One definition, used both per
template bin (makeDatacard.py) and on the full-sample sums (make_process_norms.py), so that the
normalization factors remove exactly the rate part of the template variations.

`values` is an array whose first axis runs over the members (one entry per LHEPdfWeight member):
the sum of weights (or a histogram) obtained with each member's weight.

The PDF set is identified by the number of members (UL NanoAOD):
    103  NNPDF31_nnlo_hessian_pdfas (LHA 306000-306102): 0 central, 1-100 Hessian eigenvectors,
         101 / 102 alpha_s = 0.116 / 0.120          (5FS samples: ttH, ttZ, TTTo, ttX)
    101  NNPDF31_nnlo_as_0118_nf_4 (LHA 320900-321000): 0 central, 1-100 replicas
                                                     (4FS TTbb)
Check the set of a sample with the LHEPdfWeight branch title, e.g.
    uproot.open(file)["Events"]["LHEPdfWeight"].title
"""
import numpy as np

PDF_SETS = {103: "hessian", 101: "replicas"}


def pdf_set(n_members):
    if n_members not in PDF_SETS:
        raise ValueError(f"{n_members} LHEPdfWeight members: unknown PDF set (known: {PDF_SETS})")
    return PDF_SETS[n_members]


def _relative(err, central):
    return np.nan_to_num(err / central, nan=0.0, posinf=0.0, neginf=0.0)


def pdf_up_down(values):
    """(up, down) = central member x (1 +- relative PDF uncertainty).
    hessian : sqrt(sum over the 100 eigenvectors of (member - central)^2)
    replicas: half the 16-84% interval of the 100 replicas, relative to its midpoint"""
    values = np.asarray(values, dtype=float)
    central = values[0]
    if pdf_set(len(values)) == "hessian":
        rel = _relative(np.sqrt(np.sum((values[1:101] - central) ** 2, axis=0)), central)
    else:
        q16, q84 = np.quantile(values[1:101], 0.16, axis=0), np.quantile(values[1:101], 0.84, axis=0)
        rel = _relative((q84 - q16) / 2, (q84 + q16) / 2)
    return central * (1 + rel), central * (1 - rel)


def alphas_up_down(values):
    """(up, down) = central member x (1 +- relative alpha_s uncertainty), with
    relative uncertainty = (member 102 [alpha_s 0.120] - member 101 [0.116]) / 2 / central.
    Only for the Hessian set (the 4FS replica set has no alpha_s members)."""
    values = np.asarray(values, dtype=float)
    if pdf_set(len(values)) != "hessian":
        raise ValueError("alpha_s variations exist only for the 103-member Hessian set")
    central = values[0]
    rel = _relative((values[102] - values[101]) / 2, central)
    return central * (1 + rel), central * (1 - rel)


def has_alphas(n_members):
    return pdf_set(n_members) == "hessian"
