from numba import njit


@njit
def get_partons_provenance_ttHbb_fullhadronic(pdgIds, array_builder):
    """
    Fully hadronic ttH provenance.
    Last 2 particles are always H->bb daughters appended from GenPart.

    Observed structures (N = total including appended H->bb):
    N=5:  [rad, b, -b,          H_b, H_b]  -> both tops leptonic
    N=6:  [b, -b, q, q,         H_b, H_b]  -> one top hadronic
    N=7:  [rad, b, -b, q, q,    H_b, H_b]  -> one top hadronic + radiation
    N=8:  [b, -b, q, q, q, q,   H_b, H_b]  -> BOTH tops hadronic
    N=9:  [rad, b,-b, q,q, q,q, H_b, H_b]  -> both tops hadronic + radiation

    Provenance: 1=H bb, 2=t1 b, 3=t2 b, 4=t1 W daughters, 5=t2 W daughters
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # Last 2 always H->bb appended from GenPart
        from_part[n - 2] = 1
        from_part[n - 1] = 1

        n_lhe = n - 2

        if n_lhe == 3:
            # [rad, b, -b] — both tops leptonic
            from_part[0] = -1
            if abs(ids[1]) == 5 and abs(ids[2]) == 5:
                from_part[1] = 2
                from_part[2] = 3

        elif n_lhe == 4:
            # [b, -b, q, q] — one top hadronic
            if abs(ids[0]) == 5 and abs(ids[1]) == 5:
                from_part[0] = 2
                from_part[1] = 3
                from_part[2] = 4
                from_part[3] = 4

        elif n_lhe == 5:
            # [rad, b, -b, q, q] — one top hadronic + radiation
            from_part[0] = -1
            if abs(ids[1]) == 5 and abs(ids[2]) == 5:
                from_part[1] = 2
                from_part[2] = 3
                from_part[3] = 4
                from_part[4] = 4

        elif n_lhe == 6:
            # [b, -b, q, q, q, q] — BOTH tops hadronic
            if abs(ids[0]) == 5 and abs(ids[1]) == 5:
                from_part[0] = 2
                from_part[1] = 3
                from_part[2] = 4
                from_part[3] = 4
                from_part[4] = 5
                from_part[5] = 5

        elif n_lhe == 7:
            # [rad, b, -b, q, q, q, q] — both tops hadronic + radiation
            from_part[0] = -1
            if abs(ids[1]) == 5 and abs(ids[2]) == 5:
                from_part[1] = 2
                from_part[2] = 3
                from_part[3] = 4
                from_part[4] = 4
                from_part[5] = 5
                from_part[6] = 5

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_TTZToBB_fullhadronic(pdgIds, array_builder):
    """
    Fully hadronic TTZToBB provenance.
    Last 2 particles are always Z->bb daughters appended from GenPart.

    Observed structures (N = total including appended Z->bb):
    N=7:  [b,-b, b,-b, rad,         Z_b, Z_b]  -> no W daughters
    N=8:  [b,-b, b,-b, q,q,         Z_b, Z_b]  -> one top W daughters
    N=9:  [b,-b, b,-b, q,q, rad,    Z_b, Z_b]  -> one top W + radiation
    N=10: [b,-b, b,-b, q,q, q,q,    Z_b, Z_b]  -> BOTH tops hadronic

    Structure:
    [0] = t1 b (prov=2)
    [1] = t2 b (prov=3)
    [2] = LHE Z b1 (duplicate of appended, prov=-1)
    [3] = LHE Z b2 (duplicate of appended, prov=-1)
    [4+] = W daughters and/or radiation
    [n-2,n-1] = GenPart Z->bb appended (prov=1)

    Provenance: 1=Z bb, 2=t1 b, 3=t2 b, 4=t1 W daughters, 5=t2 W daughters
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # Last 2 always Z->bb appended from GenPart
        from_part[n - 2] = 1
        from_part[n - 1] = 1

        # First 2 always top b-quarks
        if abs(ids[0]) == 5 and abs(ids[1]) == 5:
            from_part[0] = 2   # t1 b
            from_part[1] = 3   # t2 b
        # [2] and [3] are LHE Z daughters — leave as -1 (duplicates)

        n_lhe = n - 2
        # W daughters start at position 4
        n_w = n_lhe - 4  # number of partons after [0,1,2,3]

        if n_w == 1:
            # N=7: radiation only, no W daughters
            from_part[4] = -1

        elif n_w == 2:
            # N=8: one top's W daughters
            from_part[4] = 4
            from_part[5] = 4

        elif n_w == 3:
            # N=9: one top W + radiation
            from_part[4] = 4
            from_part[5] = 4
            from_part[6] = -1

        elif n_w == 4:
            # N=10: BOTH tops fully hadronic
            from_part[4] = 4
            from_part[5] = 4
            from_part[6] = 5
            from_part[7] = 5

        elif n_w == 5:
            # N=11: both tops + radiation
            from_part[4] = 4
            from_part[5] = 4
            from_part[6] = 5
            from_part[7] = 5
            from_part[8] = -1

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_ttbb4F_fullhadronic(pdgIds, array_builder):
    """
    Fully hadronic TTbb_Hadronic provenance (4F scheme).
    No H/Z appended — pure LHE partons.

    Observed structures:
    N=8:  [b,-b, b,-b, q,q, q,q]            -> g→bb + both tops hadronic (rare)
    N=9:  [b,-b, rad, b,-b, q,q, q,q]       -> g→bb + radiation + both tops hadronic (99.9%)

    Structure for N=9 (dominant):
    [0] b    → g→bb quark 1 (prov=1)
    [1] -b   → g→bb quark 2 (prov=1)
    [2] rad  → radiation (prov=-1)
    [3] b    → t1 b-quark (prov=2)
    [4] -b   → t2 b-quark (prov=3)
    [5] q    → t1 W daughter 1 (prov=4)
    [6] q    → t1 W daughter 2 (prov=4)
    [7] q    → t2 W daughter 1 (prov=5)
    [8] q    → t2 W daughter 2 (prov=5)

    Provenance: 1=g->bb, 2=t1 b, 3=t2 b, 4=t1 W daughters, 5=t2 W daughters
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # First 2 always g->bb
        from_part[0] = 1
        from_part[1] = 1

        if n == 8:
            # [b,-b, b,-b, q,q, q,q] — no radiation
            if abs(ids[2]) == 5 and abs(ids[3]) == 5:
                from_part[2] = 2   # t1 b
                from_part[3] = 3   # t2 b
                from_part[4] = 4   # t1 W daughter 1
                from_part[5] = 4   # t1 W daughter 2
                from_part[6] = 5   # t2 W daughter 1
                from_part[7] = 5   # t2 W daughter 2

        elif n == 9:
            # [b,-b, rad, b,-b, q,q, q,q] — with radiation (99.9% of events)
            from_part[2] = -1   # radiation
            if abs(ids[3]) == 5 and abs(ids[4]) == 5:
                from_part[3] = 2   # t1 b
                from_part[4] = 3   # t2 b
                from_part[5] = 4   # t1 W daughter 1
                from_part[6] = 4   # t1 W daughter 2
                from_part[7] = 5   # t2 W daughter 1
                from_part[8] = 5   # t2 W daughter 2

        # else: unexpected length, leave as -1 except g->bb

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_tt5F_fullhadronic(pdgIds, array_builder):
    """
    Fully hadronic TTToHadronic provenance (5F scheme).
    No H/Z, no g->bb.

    Observed structures:
    N=6:  [b,-b, q,q, q,q]          -> both tops hadronic, no radiation (rare, 0.001%)
    N=7:  [rad, b,-b, q,q, q,q]     -> both tops hadronic + radiation (99.9%)

    Structure for N=7 (dominant):
    [0] rad  → radiation (prov=-1)
    [1] b    → t1 b-quark (prov=2)
    [2] -b   → t2 b-quark (prov=3)
    [3] q    → t1 W daughter 1 (prov=4)
    [4] q    → t1 W daughter 2 (prov=4)
    [5] q    → t2 W daughter 1 (prov=5)
    [6] q    → t2 W daughter 2 (prov=5)

    Provenance: 2=t1 b, 3=t2 b, 4=t1 W daughters, 5=t2 W daughters
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        if n == 6:
            # [b,-b, q,q, q,q] — no radiation
            if abs(ids[0]) == 5 and abs(ids[1]) == 5:
                from_part[0] = 2   # t1 b
                from_part[1] = 3   # t2 b
                from_part[2] = 4   # t1 W daughter 1
                from_part[3] = 4   # t1 W daughter 2
                from_part[4] = 5   # t2 W daughter 1
                from_part[5] = 5   # t2 W daughter 2

        elif n == 7:
            # [rad, b,-b, q,q, q,q] — with radiation (99.9%)
            from_part[0] = -1   # radiation
            if abs(ids[1]) == 5 and abs(ids[2]) == 5:
                from_part[1] = 2   # t1 b
                from_part[2] = 3   # t2 b
                from_part[3] = 4   # t1 W daughter 1
                from_part[4] = 4   # t1 W daughter 2
                from_part[5] = 5   # t2 W daughter 1
                from_part[6] = 5   # t2 W daughter 2

        # else: unexpected length, leave as -1

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder
