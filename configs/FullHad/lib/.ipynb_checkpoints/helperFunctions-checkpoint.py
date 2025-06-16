##########################################
### Helper Functions for ZHbb analysis ###
##########################################

from numba import njit, jit, prange
import numpy as np
import pandas as pd
import awkward as ak
import math

def ak_crosscleaned(eta1, phi1, eta2, phi2, cut): # cross cleaning for arbitrary length array inputs
    # input 4 awkward arrays
    # jet is second object
    c1_ = np.array(eta1.counts)
    c2_ = np.array(eta2.counts)
    out_ = np.ones(len(eta2.flatten()))
    args_ = ak.flatten(eta1), ak.flatten(phi1), c1_, ak.flatten(eta2), ak.flatten(phi2), c2_, cut, math.pi, out_
    @njit(parallel=False)
    def njit_cc(e1,p1,c1,e2,p2,c2,cut,pi,o):
        iter1 = 0
        iter2 = 0
        a = c1.size
        for i in range(a):
            for ij in range(c2[i]):
                for ii in range(c1[i]):
                    deta = e1[ii+iter1] - e2[ij+iter2]
                    dphi = p1[ii+iter1] - p2[ij+iter2]
                    if (dphi > pi):
                        dphi = dphi - 2*pi
                    elif (dphi <= -pi):
                        dphi = dphi + 2*pi
                    dr = float((deta**2 + dphi**2)**(.5))
                    if dr <= cut:
                        o[ij+iter2] = 0
                        break
            iter1 += c1[i]
            iter2 += c2[i]
        return o
    out = ak.unflatten(njit_cc(*args_), ak.count(eta2)).astype(bool)
    return out
