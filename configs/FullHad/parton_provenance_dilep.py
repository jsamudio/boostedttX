from numba import njit


@njit
def get_partons_provenance_ttHbb_dileptonic(pdgIds, array_builder):
    """
    structures (N = total including appended H->bb):
    N=4:  [b, -b,        H_b, H_b]  -> no radiation
    N=5:  [rad, b, -b,   H_b, H_b]  -> with radiation

    Structure:
    [0]       = t1 b-quark      (ids==+5, prov=2)   [N=4]
    [1]       = t2 b-quark      (ids==-5, prov=3)   [N=4]
    [0]       = radiation       (prov=-1)            [N=5]
    [1]       = t1 b-quark      (ids==+5, prov=2)   [N=5]
    [2]       = t2 b-quark      (ids==-5, prov=3)   [N=5]
    [n-2,n-1] = H->bb appended from GenPart          (prov=1)

    Provenance:
    1  = Higgs bb (appended from GenPart)
    2  = t1 b-quark (from top)
    3  = t2 b-quark (from antitop)
    -1 = radiation
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # Last 2 always H->bb appended from GenPart
        from_part[n - 2] = 1
        from_part[n - 1] = 1

        n_lhe = n - 2  # number of LHE partons (excluding appended H->bb)

        if n_lhe == 2:
            # N=4: [b, -b] — no radiation
            if ids[0] == 5 and ids[1] == -5:
                from_part[0] = 2   # t1 b (from top)
                from_part[1] = 3   # t2 b (from antitop)

        elif n_lhe == 3:
            # N=5: [rad, b, -b] — with radiation
            from_part[0] = -1   # radiation
            if ids[1] == 5 and ids[2] == -5:
                from_part[1] = 2   # t1 b (from top)
                from_part[2] = 3   # t2 b (from antitop)

        # else: unexpected length — leave as -1 except H->bb already assigned above

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_TTZToBB_dileptonic(pdgIds, array_builder):
    """
    structures (N = total including appended Z->bb):
    N=6:  [b, -b, Z_b, Z_b,        Z_b, Z_b]  -> no radiation
    N=7:  [b, -b, Z_b, Z_b, rad,   Z_b, Z_b]  -> with radiation

    Structure:
    [0]       = t1 b-quark         (ids==+5, prov=2)
    [1]       = t2 b-quark         (ids==-5, prov=3)
    [2]       = LHE Z daughter 1   (duplicate of appended, prov=-1)
    [3]       = LHE Z daughter 2   (duplicate of appended, prov=-1)
    [4]       = radiation          (only in N=7,            prov=-1)
    [n-2,n-1] = Z->bb appended from GenPart                (prov=1)

    Provenance:
    1  = Z bb (appended from GenPart)
    2  = t1 b-quark (from top)
    3  = t2 b-quark (from antitop)
    -1 = LHE Z daughters (duplicates) or radiation
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # Last 2 always Z->bb appended from GenPart
        from_part[n - 2] = 1
        from_part[n - 1] = 1

        # [2] and [3] are LHE Z daughters — always leave as -1 (duplicates of appended)

        # First 2 always top b-quarks — check sign to distinguish top from antitop
        if ids[0] == 5 and ids[1] == -5:
            from_part[0] = 2   # t1 b (from top)
            from_part[1] = 3   # t2 b (from antitop)

        n_lhe = n - 2  # number of LHE partons (excluding appended Z->bb)

        if n_lhe == 4:
            # N=6: [b, -b, Z_b, Z_b] — no radiation, nothing extra to assign
            pass

        elif n_lhe == 5:
            # N=7: [b, -b, Z_b, Z_b, rad] — radiation at position [4]
            from_part[4] = -1   # radiation, already -1 but explicit for clarity

        # else: unexpected length — leave as -1 except assignments above

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_ttbb4F_dileptonic(pdgIds, array_builder):
    """
    Dileptonic ttbb (4-flavor scheme) provenance.
    structures:
    N=4:  [b, -b, b, -b]           -> g->bb + top b-quarks, no radiation
    N=5:  [b, -b, rad, b, -b]      -> g->bb + radiation + top b-quarks

    Structure for N=4:
    [0] b    -> g->bb quark 1  (prov=1)
    [1] -b   -> g->bb quark 2  (prov=1)
    [2] b    -> t1 b-quark     (ids==+5, prov=2)
    [3] -b   -> t2 b-quark     (ids==-5, prov=3)

    Structure for N=5:
    [0] b    -> g->bb quark 1  (prov=1)
    [1] -b   -> g->bb quark 2  (prov=1)
    [2] rad  -> radiation      (prov=-1)
    [3] b    -> t1 b-quark     (ids==+5, prov=2)
    [4] -b   -> t2 b-quark     (ids==-5, prov=3)

    Provenance:
    1  = g->bb quarks
    2  = t1 b-quark (from top)
    3  = t2 b-quark (from antitop)
    -1 = radiation
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        # First 2 always g->bb
        from_part[0] = 1
        from_part[1] = 1

        if n == 4:
            # [b, -b, b, -b] — no radiation
            if ids[2] == 5 and ids[3] == -5:
                from_part[2] = 2   # t1 b (from top)
                from_part[3] = 3   # t2 b (from antitop)

        elif n == 5:
            # [b, -b, rad, b, -b] — with radiation
            from_part[2] = -1   # radiation
            if ids[3] == 5 and ids[4] == -5:
                from_part[3] = 2   # t1 b (from top)
                from_part[4] = 3   # t2 b (from antitop)

        # else: unexpected length — leave as -1 except g->bb already assigned above

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder


@njit
def get_partons_provenance_tt5F_dileptonic(pdgIds, array_builder):
    """
    Dileptonic tt (5-flavor scheme) provenance.
    Both W bosons decay leptonically — no W decay quarks in the LHE parton array.
    No H/Z/g->bb — pure LHE partons only.

    Observed structures:
    N=2:  [b, -b]           -> no radiation (rare)
    N=3:  [rad, b, -b]      -> with radiation (dominant)

    Structure for N=2:
    [0] b    -> t1 b-quark  (ids==+5, prov=2)
    [1] -b   -> t2 b-quark  (ids==-5, prov=3)

    Structure for N=3:
    [0] rad  -> radiation   (prov=-1)
    [1] b    -> t1 b-quark  (ids==+5, prov=2)
    [2] -b   -> t2 b-quark  (ids==-5, prov=3)

    Provenance:
    2  = t1 b-quark (from top)
    3  = t2 b-quark (from antitop)
    -1 = radiation
    """
    for ids in pdgIds:
        n = len(ids)
        from_part = [-1] * n

        if n == 2:
            # [b, -b] — no radiation
            if ids[0] == 5 and ids[1] == -5:
                from_part[0] = 2   # t1 b (from top)
                from_part[1] = 3   # t2 b (from antitop)

        elif n == 3:
            # [rad, b, -b] — with radiation
            from_part[0] = -1   # radiation
            if ids[1] == 5 and ids[2] == -5:
                from_part[1] = 2   # t1 b (from top)
                from_part[2] = 3   # t2 b (from antitop)

        # else: unexpected length — leave everything as -1

        array_builder.begin_list()
        for fp in from_part:
            array_builder.append(fp)
        array_builder.end_list()

    return array_builder
