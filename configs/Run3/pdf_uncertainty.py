"""
PDF and alpha_s uncertainties from the LHEPdfWeight members. One definition, used both per
template bin (makeDatacard.py) and on the full-sample sums (make_process_norms.py), so that the
normalization factors remove exactly the rate part of the template variations.

`values` is an array whose first axis runs over the members (one entry per LHEPdfWeight member):
the sum of weights (or a histogram) obtained with each member's weight.

Run 3 version of Run2/pdf_uncertainty.py. The PDF set is identified by the number of members; the
Summer24 sets (check_lhe_weights.py, LHEPdfWeight titles, 2026-11-27) are both Hessian:
    103  NNPDF31_nnlo_as_0118_mc_hessian_pdfas (LHA 325300-325402): 0 central, 1-100 Hessian eigenvectors,
         101 / 102 alpha_s = 0.116 / 0.120     (5FS: ttH, ttZ, TTTo, TTW, TTTT, tW, THW, V+jets, QCD, SMEFT)
    101  NNPDF31_nnlo_as_0118_nf_4_mc_hessian (LHA 325500-325600): 0 central, 1-100 Hessian eigenvectors,
         no alpha_s members                    (4FS: TTbb, THQ, tZq, single top s / t channel)
Run 2 differs for the 101-member set: the UL 4FS TTbb used the replica set NNPDF31_nnlo_as_0118_nf_4
(320900), with the 16-84% replica interval. Another member count raises in pdf_set.
Check the set of a sample with check_lhe_weights.py.
"""
import numpy as np

PDF_SETS = {103: "hessian_pdfas", 101: "hessian"}


def pdf_set(n_members):
    if n_members not in PDF_SETS:
        raise ValueError(f"{n_members} LHEPdfWeight members: unknown PDF set (known: {PDF_SETS})")
    return PDF_SETS[n_members]


def _relative(err, central):
    return np.nan_to_num(err / central, nan=0.0, posinf=0.0, neginf=0.0)


def pdf_up_down(values):
    """(up, down) = central member x (1 +- relative PDF uncertainty): sqrt of the sum over the 100
    Hessian eigenvectors (members 1-100) of (member - central)^2, for both Summer24 sets."""
    values = np.asarray(values, dtype=float)
    pdf_set(len(values))
    central = values[0]
    rel = _relative(np.sqrt(np.sum((values[1:101] - central) ** 2, axis=0)), central)
    return central * (1 + rel), central * (1 - rel)


def alphas_up_down(values):
    """(up, down) = central member x (1 +- relative alpha_s uncertainty), with
    relative uncertainty = (member 102 [alpha_s 0.120] - member 101 [0.116]) / 2 / central.
    Only for the 103-member set (the 4FS set has no alpha_s members)."""
    values = np.asarray(values, dtype=float)
    if pdf_set(len(values)) != "hessian_pdfas":
        raise ValueError("alpha_s variations exist only for the 103-member set")
    central = values[0]
    rel = _relative((values[102] - values[101]) / 2, central)
    return central * (1 + rel), central * (1 - rel)


def has_alphas(n_members):
    return pdf_set(n_members) == "hessian_pdfas"
