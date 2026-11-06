#!/usr/bin/env python3
"""
Generator-level check of the DCTR hdamp weights on the 4FS ttbb (no selection, NanoAOD directly).

DCTR is trained on the gen top / antitop kinematics, so with the right convention the reweighted
nominal sample must reproduce the gen top kinematics of the dedicated hdamp samples. This tests
(a) the weight convention, p0/p1 (as in configs/common/dctr.py) or p1/p0, and (b) how well the
ttbar-trained classifier carries over to ttbb, with much more statistics than the analysis selection.

Per sample (nominal, hdampUP, hdampDOWN) the first --max-events events of the dataset JSON files
are read (GenPart first-copy t / tbar, sign(genWeight)). The features and the model call are those
of dctr.py. Histograms (unit area) of pT(tt), m(tt), y(tt), top pT; ratio to nominal for the sample,
DCTR p0/p1 and DCTR p1/p0; chi2 of each DCTR convention against the sample (stat. of both).
Also the inclusive mean weight (1 expected for a normalization-preserving reweighting).

    python3 dctr_gen_check.py --sample TTbb_SemiLeptonic --era 2018 --max-events 300000
"""
import argparse
import json
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
import dctr  # noqa: E402
from inference_common import DCTR_HDAMP  # noqa: E402

FIRST_COPY_BIT = 12          # GenPart_statusFlags bit isFirstCopy
OBS = {
    "pt_tt": (np.array([0, 25, 50, 75, 100, 150, 200, 300, 400, 600, 1000]), r"$p_T(t\bar{t})$ [GeV]"),
    "m_tt": (np.array([340, 400, 450, 500, 600, 700, 850, 1000, 1500, 3000]), r"$m(t\bar{t})$ [GeV]"),
    "y_tt": (np.linspace(-2.5, 2.5, 11), r"$y(t\bar{t})$"),
    "pt_top": (np.array([0, 50, 100, 150, 200, 300, 400, 600, 1000]), r"top $p_T$ [GeV]"),
}


def files_of(name, era):
    d = json.load(open(f"datasets/{name}.json"))
    key = f"{name}_{era}"
    return [f.replace("file:", "") for f in d[key]["files"]]


def read(files, max_events):
    """Arrays of t and tbar (pt, eta, phi, mass) for events with exactly one first-copy t and tbar."""
    cols = ["GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_pdgId",
            "GenPart_statusFlags", "genWeight"]
    out, n = [], 0
    for path in files:
        if n >= max_events:
            break
        a = uproot.open(path)["Events"].arrays(cols, entry_stop=max_events - n)
        n += len(a)
        first = (a.GenPart_statusFlags >> FIRST_COPY_BIT) & 1 == 1
        sel = {}
        for name, pid in (("t", 6), ("tbar", -6)):
            m = first & (a.GenPart_pdgId == pid)
            sel[name] = m
        ok = (ak.sum(sel["t"], axis=1) == 1) & (ak.sum(sel["tbar"], axis=1) == 1)
        ev = {"w": np.sign(ak.to_numpy(a.genWeight[ok]))}
        for name in ("t", "tbar"):
            for v in ("pt", "eta", "phi", "mass"):
                ev[f"{name}_{v}"] = ak.to_numpy(ak.flatten(a[f"GenPart_{v}"][sel[name]][ok]))
        out.append(ev)
        print(f"  {os.path.basename(path)}: {len(a)} events, {int(ak.sum(ok))} with one t and one tbar")
    return {k: np.concatenate([o[k] for o in out]) for k in out[0]}


def p4(ev, name):
    pt, eta, phi, m = (ev[f"{name}_{v}"] for v in ("pt", "eta", "phi", "mass"))
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return px, py, pz, np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2)


def rapidity(pz, e):
    return 0.5 * np.log((e + pz) / (e - pz))


def observables(ev):
    t, tb = p4(ev, "t"), p4(ev, "tbar")
    px, py, pz, e = (a + b for a, b in zip(t, tb))
    return {"pt_tt": np.hypot(px, py), "m_tt": np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0)),
            "y_tt": rapidity(pz, e), "pt_top": ev["t_pt"]}


# hdamp / mt of the UL CP5 samples: nominal 1.379, up 2.305, down 0.8738 (hdamp = 1.379 +0.926 -0.505 mt)
HDAMP_TARGET = {"up": 2.305, "down": 0.8738}


def features(ev, name, pid, hdamp):
    px, py, pz, e = p4(ev, name)
    n = len(px)
    return np.stack([np.log10(ev[f"{name}_pt"]), rapidity(pz, e), ev[f"{name}_phi"],
                     ev[f"{name}_mass"] / dctr.MAX_MASS, np.full(n, pid), np.full(n, hdamp)], axis=-1)


def dctr_ratio(ev, path, pt_tt, hdamp=dctr.HDAMP_NOMINAL):
    """p0/p1 per event as in dctr.py (1 above MAX_TTBAR_PT); hdamp = value of the 6th feature."""
    inputs = np.stack([features(ev, "t", 0.1, hdamp), features(ev, "tbar", 0.2, hdamp)], axis=1).astype(np.float32)
    r = np.ones(len(pt_tt))
    for i in range(0, len(r), 200000):
        r[i:i + 200000] = dctr._ratio(path, inputs[i:i + 200000])
    return np.where(pt_tt < dctr.MAX_TTBAR_PT, r, 1.0)


def hist(x, w, bins):
    x = np.clip(x, bins[0], bins[-1] - 1e-9)
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w ** 2)
    s = h.sum()
    return h / s, np.sqrt(e2) / s


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", default="TTbb_SemiLeptonic")
    ap.add_argument("--era", default="2018")
    ap.add_argument("--max-events", type=int, default=300000)
    ap.add_argument("--outdir", default="plots/hdamp_dctr_gen")
    ap.add_argument("--hdamp-feature", choices=["nominal", "target"], default="nominal",
                    help="6th DCTR feature: nominal hdamp/mt 1.379 for both models (as dctr.py; the TOP MLReweighting "
                         "twiki convention), or the target value of each model (up 2.305, down 0.8738; NOT the "
                         "twiki convention, kept only as a cross-check)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    store = {"sample": args.sample, "era": args.era}       # twiki convention (p0/p1), for plot_dctr_summary.py

    evs = {}
    for tag, name in (("nominal", args.sample), ("up", f"{args.sample}_hdampUP"), ("down", f"{args.sample}_hdampDOWN")):
        print(f"Reading {name} {args.era}")
        evs[tag] = read(files_of(name, args.era), args.max_events)
    obs = {k: observables(v) for k, v in evs.items()}
    nom, w0 = evs["nominal"], evs["nominal"]["w"]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for var in ("up", "down"):
        hd_feat = dctr.HDAMP_NOMINAL if args.hdamp_feature == "nominal" else HDAMP_TARGET[var]
        r = dctr_ratio(nom, DCTR_HDAMP["models"][var], obs["nominal"]["pt_tt"], hd_feat)
        conv = {"p0/p1": r, "p1/p0": 1.0 / r}
        print(f"\n=== hdamp {var}: {len(nom['w'])} nominal, {len(evs[var]['w'])} {var} events, "
              f"hdamp feature {hd_feat}")
        store[f"{var}_meanw"] = np.sum(w0 * r) / np.sum(w0)
        for c, w in conv.items():
            print(f"  DCTR {c}: inclusive mean weight {np.sum(w0 * w) / np.sum(w0):.4f}, "
                  f"min {w.min():.3f} max {w.max():.3f}")
        for o, (bins, label) in OBS.items():
            hn, en = hist(obs["nominal"][o], w0, bins)
            if o == "pt_tt":      # direction of the effect: sample and DCTR ratio per bin
                hs_, _ = hist(obs[var][o], evs[var]["w"], bins)
                print("    pt_tt ratio to nominal per bin  sample: " + " ".join(f"{v:.3f}" for v in hs_ / hn))
                for c, w in conv.items():
                    hd_, _ = hist(obs["nominal"][o], w0 * w, bins)
                    print(f"    pt_tt ratio to nominal per bin  {c:6}: " + " ".join(f"{v:.3f}" for v in hd_ / hn))
            hs, es = hist(obs[var][o], evs[var]["w"], bins)
            ok = hn > 0
            rs = np.divide(hs, hn, out=np.zeros_like(hn), where=ok)
            line, fig = [], plt.figure(figsize=(6, 4))
            ax = fig.gca()
            x = 0.5 * (bins[1:] + bins[:-1])
            ax.errorbar(x, rs, yerr=np.divide(es, hn, out=np.zeros_like(hn), where=ok), fmt="ko", ms=4,
                        label=f"hdamp{var.upper()} sample / nominal")
            for c, w in conv.items():
                hd, ed = hist(obs["nominal"][o], w0 * w, bins)
                chi2 = (((hs - hd) ** 2) / np.maximum(es ** 2 + ed ** 2, 1e-30))[ok & (hs > 0)].sum()
                line.append(f"{c} {chi2:.1f}/{int((ok & (hs > 0)).sum()) - 1}")
                if c == "p0/p1":
                    store.update({f"{var}_{o}_bins": bins, f"{var}_{o}_rs": rs,
                                  f"{var}_{o}_es": np.divide(es, hn, out=np.zeros_like(hn), where=ok),
                                  f"{var}_{o}_rd": np.divide(hd, hn, out=np.zeros_like(hn), where=ok),
                                  f"{var}_{o}_ed": np.divide(ed, hn, out=np.zeros_like(hn), where=ok),
                                  f"{var}_{o}_chi2": chi2, f"{var}_{o}_ndf": int((ok & (hs > 0)).sum()) - 1})
                ax.stairs(np.divide(hd, hn, out=np.zeros_like(hn), where=ok), bins, baseline=None, lw=1.5,
                          label=f"DCTR {c} / nominal")
            print(f"    {o:8} chi2/ndf (sample vs DCTR): " + "  ".join(line))
            ax.axhline(1, color="gray", lw=0.8)
            ax.set_xlabel(label)
            ax.set_ylabel("ratio to nominal (unit area)")
            ax.legend(fontsize=7)
            ax.set_title(f"{args.sample} {args.era}, gen level: hdamp {var}", fontsize=9)
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(f"{args.outdir}/{args.sample}_{args.era}_{var}_{o}_hd{args.hdamp_feature}.{ext}")
            plt.close(fig)

    if args.hdamp_feature == "nominal":
        np.savez(f"{args.outdir}/{args.sample}_{args.era}_dctr.npz", **store)
        print(f"Wrote {args.outdir}/{args.sample}_{args.era}_dctr.npz (for plot_dctr_summary.py)")
