#!/usr/bin/env python3
"""
Generator-level check of the DCTR hdamp weights on the 4FS ttbb (no selection, NanoAOD directly).

Also the additional b system (b jets not from the top decays), the part of the event that hdamp does not
act on in the 4FS ttbb (bb in the Born) but does in the 5FS (bb from the hardest emission / shower).

DCTR is trained on the gen top / antitop kinematics, so with the right convention the reweighted
nominal sample must reproduce the gen top kinematics of the dedicated hdamp samples. This tests
(a) the weight convention, p0/p1 (as in configs/common/dctr.py) or p1/p0, and (b) how well the
ttbar-trained classifier carries over to ttbb, with much more statistics than the analysis selection.

Per sample (nominal, hdampUP, hdampDOWN) the first --max-events events of the dataset JSON files
are read (GenPart first-copy t / tbar, sign(genWeight)). The features and the model call are those
of dctr.py. Histograms (unit area) of pT(tt), m(tt), y(tt), top pT; ratio to nominal for the sample,
DCTR p0/p1 and DCTR p1/p0; chi2 of each DCTR convention against the sample (stat. of both).
Additional b jets: GenJet pT > 20 GeV, |eta| < 2.4, hadronFlavour 5, dR > 0.4 to both b quarks from the
top decays (GenPart |pdgId| 5 with a top mother), as in the genTtbarId tt+B definition. Observables:
number of additional b jets (all events); pT(bb), m(bb), dR(bb), pT(ttbb) (two leading additional b
jets + t + tbar: the recoil against the extra radiation, which hdamp controls in the 4FS) for events
with >= 2; leading additional b-jet pT for events with >= 1. Each histogram is unit area over the
events where the observable is defined.
Additional light jets: GenJet pT > 20 GeV, |eta| < 2.4, hadronFlavour != 5 (udscg), dR > 0.4 to every
top-decay product (b quarks with a top mother, quarks and charged leptons with a W mother; GenJets
include leptons). Observables: their number (all events) and the leading one's pT (events with >= 1),
the extra radiation that hdamp acts on directly in both schemes.
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
ADDB_PT, ADDB_ETA, ADDB_DR = 20.0, 2.4, 0.4
OBS = {
    "pt_tt": (np.array([0, 25, 50, 75, 100, 150, 200, 300, 400, 600, 1000]), r"$p_T(t\bar{t})$ [GeV]"),
    "m_tt": (np.array([340, 400, 450, 500, 600, 700, 850, 1000, 1500, 3000]), r"$m(t\bar{t})$ [GeV]"),
    "y_tt": (np.linspace(-2.5, 2.5, 11), r"$y(t\bar{t})$"),
    "pt_top": (np.array([0, 50, 100, 150, 200, 300, 400, 600, 1000]), r"top $p_T$ [GeV]"),
    "n_addb": (np.arange(-0.5, 5.5, 1.0), "N additional b jets"),
    "pt_b1": (np.array([20, 30, 40, 60, 80, 100, 150, 200, 300, 500]), r"leading additional b jet $p_T$ [GeV]"),
    "pt_bb": (np.array([0, 25, 50, 75, 100, 150, 200, 300, 500, 1000]), r"$p_T(b\bar{b})$, additional [GeV]"),
    "m_bb": (np.array([0, 30, 60, 90, 120, 160, 200, 300, 500, 1000]), r"$m(b\bar{b})$, additional [GeV]"),
    "dr_bb": (np.linspace(0.4, 4.4, 9), r"$\Delta R(b\bar{b})$, additional"),
    "n_addj": (np.arange(-0.5, 6.5, 1.0), "N additional light jets"),
    "pt_j1": (np.array([20, 30, 40, 60, 80, 100, 150, 200, 300, 500, 1000]), r"leading additional light jet $p_T$ [GeV]"),
    "pt_ttbb": (np.array([0, 25, 50, 75, 100, 150, 200, 300, 400, 600, 1000]), r"$p_T(t\bar{t}b\bar{b})$ [GeV]"),
}


def files_of(name, era):
    d = json.load(open(f"datasets/{name}.json"))
    key = f"{name}_{era}"
    return [f.replace("file:", "") for f in d[key]["files"]]


def read(files, max_events, ttb=False):
    """Arrays of t and tbar (pt, eta, phi, mass) for events with exactly one first-copy t and tbar; with
    ttb=True only tt+B events (genTtbarId % 100 >= 51). Read in steps (GenPart is large). Also returns the
    sign(genWeight) sums of all events and of the tt+B events (tt+B fraction)."""
    cols = ["GenPart_pt", "GenPart_eta", "GenPart_phi", "GenPart_mass", "GenPart_pdgId",
            "GenPart_statusFlags", "GenPart_genPartIdxMother", "genWeight",
            "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_mass", "GenJet_hadronFlavour"] \
        + (["genTtbarId"] if ttb else [])
    out, n, sum_all, sum_ttb = [], 0, 0.0, 0.0
    for path in files:
        if n >= max_events:
            break
        tree = uproot.open(path)["Events"]
        n_file = 0
        for a in tree.iterate(cols, step_size=200000, entry_stop=min(tree.num_entries, max_events - n)):
            n_file += len(a)
            first = (a.GenPart_statusFlags >> FIRST_COPY_BIT) & 1 == 1
            sel = {name: first & (a.GenPart_pdgId == pid) for name, pid in (("t", 6), ("tbar", -6))}
            ok = (ak.sum(sel["t"], axis=1) == 1) & (ak.sum(sel["tbar"], axis=1) == 1)
            sign = np.sign(ak.to_numpy(a.genWeight))
            sum_all += sign.sum()
            if ttb:
                is_ttb = ak.to_numpy(a.genTtbarId % 100 >= 51)
                sum_ttb += sign[is_ttb].sum()
                ok = ok & is_ttb
            ev = {"w": sign[ak.to_numpy(ok)]}
            for name in ("t", "tbar"):
                for v in ("pt", "eta", "phi", "mass"):
                    ev[f"{name}_{v}"] = ak.to_numpy(ak.flatten(a[f"GenPart_{v}"][sel[name]][ok]))
            ev.update(additional_b(a[ok]))
            ev.update(additional_light(a[ok]))
            out.append(ev)
        n += n_file
        print(f"  {os.path.basename(path)}: {n_file} events read")
    ev = {k: np.concatenate([o[k] for o in out]) for k in out[0]}
    print(f"  -> {len(ev['w'])} events kept" + (f" (tt+B; tt+B fraction {sum_ttb / sum_all:.4f})" if ttb else ""))
    return ev, (sum_all, sum_ttb)


def additional_b(a):
    """Additional b jets per event (see the docstring): their number and the (pt, eta, phi, mass) of the
    two leading ones (NaN if absent)."""
    mom = a.GenPart_genPartIdxMother
    mom_pdg = a.GenPart_pdgId[ak.where(mom >= 0, mom, 0)]
    btop = (abs(a.GenPart_pdgId) == 5) & (mom >= 0) & (abs(mom_pdg) == 6)
    bq_eta, bq_phi = a.GenPart_eta[btop], a.GenPart_phi[btop]
    jsel = (a.GenJet_pt > ADDB_PT) & (abs(a.GenJet_eta) < ADDB_ETA) & (a.GenJet_hadronFlavour == 5)
    jets = {v: a[f"GenJet_{v}"][jsel] for v in ("pt", "eta", "phi", "mass")}
    j, b = ak.unzip(ak.cartesian([ak.zip({"eta": jets["eta"], "phi": jets["phi"]}),
                                  ak.zip({"eta": bq_eta, "phi": bq_phi})], nested=True))    # [event][jet][b]
    dphi = (j.phi - b.phi + np.pi) % (2 * np.pi) - np.pi
    dr_min = ak.fill_none(ak.min(np.hypot(j.eta - b.eta, dphi), axis=2), 99.0)
    add = dr_min > ADDB_DR
    out = {"n_addb": ak.to_numpy(ak.sum(add, axis=1)).astype(float)}
    for v in ("pt", "eta", "phi", "mass"):
        lead = ak.pad_none(jets[v][add], 2, clip=True)       # GenJets are pT-ordered
        for k in (0, 1):
            out[f"b{k + 1}_{v}"] = ak.to_numpy(ak.fill_none(lead[:, k], np.nan)).astype(float)
    return out


def min_dr(eta1, phi1, eta2, phi2):
    """Per object of collection 1, the smallest dR to collection 2 (99 if collection 2 is empty)."""
    j, b = ak.unzip(ak.cartesian([ak.zip({"eta": eta1, "phi": phi1}), ak.zip({"eta": eta2, "phi": phi2})],
                                 nested=True))
    dphi = (j.phi - b.phi + np.pi) % (2 * np.pi) - np.pi
    return ak.fill_none(ak.min(np.hypot(j.eta - b.eta, dphi), axis=2), 99.0)


def additional_light(a):
    """Additional non-b jets (see the docstring): their number and the leading one's pT (NaN if none)."""
    mom = a.GenPart_genPartIdxMother
    mom_pdg = abs(a.GenPart_pdgId[ak.where(mom >= 0, mom, 0)])
    pdg = abs(a.GenPart_pdgId)
    decay = (mom >= 0) & (((pdg == 5) & (mom_pdg == 6))
                          | ((mom_pdg == 24) & ((pdg <= 5) | (pdg == 11) | (pdg == 13) | (pdg == 15))))
    jsel = (a.GenJet_pt > ADDB_PT) & (abs(a.GenJet_eta) < ADDB_ETA) & (a.GenJet_hadronFlavour != 5)
    pt, eta, phi = a.GenJet_pt[jsel], a.GenJet_eta[jsel], a.GenJet_phi[jsel]
    add = min_dr(eta, phi, a.GenPart_eta[decay], a.GenPart_phi[decay]) > ADDB_DR
    lead = ak.pad_none(pt[add], 1, clip=True)[:, 0]        # GenJets are pT-ordered
    return {"n_addj": ak.to_numpy(ak.sum(add, axis=1)).astype(float),
            "j1_pt": ak.to_numpy(ak.fill_none(lead, np.nan)).astype(float)}


def p4(ev, name):
    pt, eta, phi, m = (ev[f"{name}_{v}"] for v in ("pt", "eta", "phi", "mass"))
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    return px, py, pz, np.sqrt(px ** 2 + py ** 2 + pz ** 2 + m ** 2)


def rapidity(pz, e):
    return 0.5 * np.log((e + pz) / (e - pz))


def observables(ev):
    t, tb = p4(ev, "t"), p4(ev, "tbar")
    px, py, pz, e = (a + b for a, b in zip(t, tb))
    out = {"pt_tt": np.hypot(px, py), "m_tt": np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0)),
           "y_tt": rapidity(pz, e), "pt_top": ev["t_pt"]}
    # additional b system: NaN (dropped from the histogram) where fewer than 2 (1 for pt_b1) b jets
    b1, b2 = p4(ev, "b1"), p4(ev, "b2")
    bx, by, bz, be = (a + b for a, b in zip(b1, b2))
    dphi = (ev["b1_phi"] - ev["b2_phi"] + np.pi) % (2 * np.pi) - np.pi
    out.update({"n_addj": ev["n_addj"], "pt_j1": ev["j1_pt"]})
    out.update({"n_addb": ev["n_addb"], "pt_b1": ev["b1_pt"], "pt_bb": np.hypot(bx, by),
                "m_bb": np.sqrt(np.maximum(be ** 2 - bx ** 2 - by ** 2 - bz ** 2, 0)),
                "dr_bb": np.hypot(ev["b1_eta"] - ev["b2_eta"], dphi),
                "pt_ttbb": np.hypot(px + bx, py + by)})
    two = ev["n_addb"] >= 2
    for k in ("pt_bb", "m_bb", "dr_bb", "pt_ttbb"):
        out[k] = np.where(two, out[k], np.nan)
    return out


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
    fin = np.isfinite(x)                       # observable undefined (e.g. < 2 additional b jets)
    x, w = np.clip(x[fin], bins[0], bins[-1] - 1e-9), w[fin]
    h, _ = np.histogram(x, bins, weights=w)
    e2, _ = np.histogram(x, bins, weights=w ** 2)
    s = h.sum()
    return h / s, np.sqrt(e2) / s


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", default="TTbb_SemiLeptonic")
    ap.add_argument("--era", default="2018")
    ap.add_argument("--max-events", type=int, default=300000)
    ap.add_argument("--ttb", action="store_true",
                    help="only tt+B events (genTtbarId %% 100 >= 51), e.g. the tt+B part of the 5FS ttbar "
                         "(~1-2%% of the events: use a few million --max-events); prints the tt+B fraction change")
    ap.add_argument("--outdir", default="plots/hdamp_dctr_gen")
    ap.add_argument("--hdamp-feature", choices=["nominal", "target"], default="nominal",
                    help="6th DCTR feature: nominal hdamp/mt 1.379 for both models (as dctr.py; the TOP MLReweighting "
                         "twiki convention), or the target value of each model (up 2.305, down 0.8738; NOT the "
                         "twiki convention, kept only as a cross-check)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sums = {}
    store = {"sample": args.sample + (" tt+B" if args.ttb else ""), "era": args.era}       # twiki convention (p0/p1), for plot_dctr_summary.py

    evs = {}
    for tag, name in (("nominal", args.sample), ("up", f"{args.sample}_hdampUP"), ("down", f"{args.sample}_hdampDOWN")):
        print(f"Reading {name} {args.era}")
        evs[tag], sums[tag] = read(files_of(name, args.era), args.max_events, args.ttb)
    if args.ttb:
        f0 = sums["nominal"][1] / sums["nominal"][0]
        for v in ("up", "down"):
            fv = sums[v][1] / sums[v][0]
            store[f"{v}_ttbfrac_ratio"] = fv / f0
            print(f"tt+B fraction: nominal {f0:.4f}, hdamp {v} {fv:.4f}, ratio {fv / f0:.4f}")
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
            if o in ("pt_tt", "pt_ttbb", "pt_bb", "pt_j1"):      # direction of the effect: sample and DCTR ratio per bin
                hs_, _ = hist(obs[var][o], evs[var]["w"], bins)
                print(f"    {o} ratio to nominal per bin  sample: " + " ".join(f"{v:.3f}" for v in hs_ / hn))
                for c, w in conv.items():
                    hd_, _ = hist(obs["nominal"][o], w0 * w, bins)
                    print(f"    {o} ratio to nominal per bin  {c:6}: " + " ".join(f"{v:.3f}" for v in hd_ / hn))
            hs, es = hist(obs[var][o], evs[var]["w"], bins)
            # unit-area histograms and number of events where the observable is defined
            # (plot_hdamp_4fs_5fs_gen.py: nominal shapes, hdamp ratios with the stat. of both samples)
            store.update({f"nominal_{o}_h": hn, f"nominal_{o}_e": en, f"{var}_{o}_h": hs, f"{var}_{o}_e": es,
                          f"nominal_{o}_n": int(np.isfinite(obs["nominal"][o]).sum()),
                          f"{var}_{o}_n": int(np.isfinite(obs[var][o]).sum())})
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
                fig.savefig(f"{args.outdir}/{args.sample}_{args.era}{'_ttB' if args.ttb else ''}_{var}_{o}_hd{args.hdamp_feature}.{ext}")
            plt.close(fig)

    if args.hdamp_feature == "nominal":
        stem = f"{args.outdir}/{args.sample}_{args.era}{'_ttB' if args.ttb else ''}_dctr.npz"
        np.savez(stem, **store)
        print(f"Wrote {stem} (for plot_dctr_summary.py)")
