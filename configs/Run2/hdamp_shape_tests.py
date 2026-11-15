"""
Shape tests of a two-sided variation (up, down) against a nominal, per datacard bin (numpy only).
Used by compare_hdamp_dctr.py (hdamp samples) and compare_ttb_theory.py (card nuisances).

All ratios are shape-only: each variation is scaled by its rate R = sum(var) / sum(nom) (one common
rate, the lnN hypothesis), so what remains is what a free normalization cannot absorb.

  r_up, r_down  (var_i / R) / nom_i                         shape ratio per bin
  A             (r_up - r_down) / 2                          hdamp-ordered (antisymmetric) part: what
                                                             moves up and down in opposite directions.
                                                             A shared nominal cancels in it, so its error
                                                             is the statistics of the two variation
                                                             samples (plus a small nominal term ~ A).
  S             (r_up + r_down) / 2 - 1                      common (symmetric) part: what moves both the
                                                             same way, e.g. a fluctuation of the shared
                                                             nominal; error includes the nominal stat.

Errors assume independent variation samples (dedicated samples like hdamp). For weight variations
of the same events (PDF, scales, PS) the per-bin values are effects, not tests: use them for sizes.
"""
import numpy as np


def _div(a, b, fill=0.0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.divide(a, b, out=np.full_like(b, fill), where=b > 0)


def collapse(h, e, n_nn):
    """Sum the mass bins of a channel flattened NN-major (bin = i_nn * n_mass + i_mass) -> per NN bin."""
    h, e = np.asarray(h, float), np.asarray(e, float)
    if len(h) % n_nn:
        raise ValueError(f"{len(h)} bins are not {n_nn} NN bins x mass bins")
    n_mass = len(h) // n_nn
    return h.reshape(n_nn, n_mass).sum(1), np.sqrt((e.reshape(n_nn, n_mass) ** 2).sum(1))


def split(hn, en, hu, eu, hd, ed, rate_u=None, rate_d=None):
    """Shape ratios and their antisymmetric / symmetric parts, with errors (independent samples).
    Rates default to the totals of the given bins."""
    hn, en, hu, eu, hd, ed = (np.asarray(x, float) for x in (hn, en, hu, eu, hd, ed))
    rate_u = hu.sum() / hn.sum() if rate_u is None else rate_u
    rate_d = hd.sum() / hn.sum() if rate_d is None else rate_d
    rn = _div(en, hn)                                   # relative nominal stat.
    ru, rd = _div(hu / rate_u, hn, 1.0), _div(hd / rate_d, hn, 1.0)
    su, sd = _div(eu / rate_u, hn), _div(ed / rate_d, hn)     # variation stat. relative to nominal
    a, s = (ru - rd) / 2, (ru + rd) / 2 - 1
    out = {"rate_up": rate_u, "rate_down": rate_d, "r_up": ru, "r_down": rd,
           "e_up": np.sqrt(su ** 2 + (ru * rn) ** 2), "e_down": np.sqrt(sd ** 2 + (rd * rn) ** 2),
           "A": a, "eA": np.sqrt((su ** 2 + sd ** 2) / 4 + (a * rn) ** 2),
           "S": s, "eS": np.sqrt((su ** 2 + sd ** 2) / 4 + ((s + 1) * rn) ** 2),
           "ok": (hn > 0) & (hu > 0) & (hd > 0)}
    return out


def chi2(x, ex, ok):
    """sum (x / ex)^2 over usable bins, and their number."""
    ok = ok & (ex > 0)
    return float(np.sum((x[ok] / ex[ok]) ** 2)), int(ok.sum())


def slope(x, y, ey, ok=None):
    """Weighted straight-line fit y = a + b (x - <x>): (b, error of b)."""
    x, y, ey = (np.asarray(v, float) for v in (x, y, ey))
    ok = (ey > 0) if ok is None else (ok & (ey > 0))
    w = 1 / ey[ok] ** 2
    xm = np.sum(w * x[ok]) / w.sum()
    sxx = np.sum(w * (x[ok] - xm) ** 2)
    return float(np.sum(w * (x[ok] - xm) * y[ok]) / sxx), float(1 / np.sqrt(sxx))


def rms(x, weights, ok=None):
    """Yield-weighted RMS of x."""
    ok = np.ones(len(x), bool) if ok is None else ok
    w = np.where(ok, np.asarray(weights, float), 0)
    return float(np.sqrt(np.sum(w * np.asarray(x) ** 2) / w.sum()))


def pvalue(c2, ndf):
    try:
        from scipy.stats import chi2 as _chi2
        return f"{_chi2.sf(c2, ndf):.2f}"
    except ImportError:
        return "n/a"
