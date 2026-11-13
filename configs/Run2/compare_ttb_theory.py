#!/usr/bin/env python3
"""
Size of every theory nuisance on one process (default tt_B) in a datacard, next to hdamp from the
dedicated 4FS samples: is hdamp_ttbb small compared with PDF, scales, parton shower, top pT, ...?

Inputs: a datacard (.txt + its .root; per year or the summed Run 2 card of sum_datacards.py) and,
optionally, the per-datacard-bin hdamp yields of compare_hdamp_dctr.py --export (same model and
years). For each nuisance acting on the process (lnN and shape, theory group by default):

  rate      Up / nominal and Down / nominal summed over all channels (lnN: its kappas), and per channel.
            Note: tt_B floats (CMS_ttbbnorm rateParam), so a common rate is absorbed by the fit.
  shape     what the free normalization cannot absorb: the variation scaled to the nominal total,
            s_i = (Up_i / rate) / Nom_i - 1 over all bins of all channels, summarized as
              rms    yield-weighted RMS of s_i (%)
              max    largest |s_i| in bins with nominal stat. < 20%
              chi2   sum of (Nom_i s_i / sigma_i)^2 with sigma_i = nominal MC stat. of the bin: whether the
                     shape is resolved beyond the template stat. (ndf = number of bins - 1)
  hdamp (4FS samples): the same, from independent samples, so their stat. fluctuations look like shape;
            'rms (stat-subtracted)' removes the expected stat. RMS of the ratio in quadrature and chi2 uses
            the stat. of both samples (chi2 ~ ndf = no shape beyond the sample stat.).

Writes a table (stdout and <outdir>/<process>_theory.csv), a bar chart of rate and shape size per
nuisance (<process>_theory_summary.png) and the per-bin shape-only ratios of the largest nuisances
(<process>_theory_bins.png).

    python3 compare_hdamp_dctr.py --no-dctr --nominal outputs/nominal/bkg/output_all.coffea \\
        --hdamp outputs_hdamp_run2/nominal/variations/output_all.coffea \\
        --edges nn_edges_disc_<model>.json --eras 2016_PreVFP 2016_PostVFP 2017 2018 --export --outdir plots/hdamp_run2
    python3 compare_ttb_theory.py datacards/datacard_nominal_disc_<model>_eftTTB_Run2.txt \\
        --hdamp plots/hdamp_run2/hdamp_ttB_<model>_2016_PostVFP_2016_PreVFP_2017_2018.json
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "common"))
from sum_datacards import parse_card, Templates, lnn  # noqa: E402

MAX_REL_STAT = 0.2          # bins used for 'max': nominal relative stat. below this


def shape_metrics(nom, enom, up, eup=None):
    """Rate and shape-only size of one variation over the concatenated bins of all channels.
    eup: stat. of an independent variation sample (hdamp); None for weight variations of the same events."""
    rate = up.sum() / nom.sum()
    ok = nom > 0
    s = np.zeros_like(nom)
    s[ok] = up[ok] / rate / nom[ok] - 1
    w = np.where(ok, nom, 0)
    rms = np.sqrt(np.sum(w * s ** 2) / w.sum())
    rel = np.divide(enom, nom, out=np.full_like(nom, np.inf), where=ok)
    good = ok & (rel < MAX_REL_STAT)
    mx = np.max(np.abs(s[good])) if good.any() else np.nan
    out = {"rate": rate, "rms": rms, "max": mx, "s": s}
    if eup is None:
        out["chi2"] = np.sum((s[ok] / rel[ok]) ** 2)
    else:
        # relative stat. of the shape ratio: both samples, the variation scaled by 1/rate
        rv = np.divide(eup, up, out=np.zeros_like(up), where=up > 0)
        er = np.sqrt(rel ** 2 + rv ** 2)
        both = ok & (up > 0)
        out["chi2"] = np.sum((s[both] / er[both]) ** 2)
        noise = np.sum(w[both] * er[both] ** 2) / w[both].sum()
        out["rms_sub"] = np.sqrt(max(rms ** 2 - noise, 0.0))
        out["rms_noise"] = np.sqrt(noise)          # RMS expected from the sample stat. alone
    out["ndf"] = int(ok.sum()) - 1
    return out


def card_nuisances(card, process, groups):
    """{name: (type, {channel: value})} of the nuisances with an entry for the process."""
    keep = None
    if groups:
        keep = {n for g in groups for n in card["groups"].get(g, [])}
    out = {}
    for name, (typ, values) in card["nuisances"].items():
        if keep is not None and name not in keep:
            continue
        entries = {ch: v for (ch, p, _), v in zip(card["columns"], values) if p == process and v != "-"}
        if entries:
            out[name] = (typ, entries)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("card", help="datacard .txt (its .root next to it)")
    ap.add_argument("--process", default="tt_B")
    ap.add_argument("--groups", nargs="*", default=["theo"],
                    help="nuisance groups to include (default: theo); empty = every nuisance on the process")
    ap.add_argument("--hdamp", default=None, help="compare_hdamp_dctr.py --export JSON (hdamp from the 4FS samples)")
    ap.add_argument("--top", type=int, default=6, help="nuisances drawn per bin (largest shape RMS)")
    ap.add_argument("--outdir", default="plots/ttb_theory")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    card = parse_card(args.card)
    tpl = Templates(card)
    channels = list(dict.fromkeys(ch for ch, p, _ in card["columns"] if p == args.process))
    nominal = {ch: tpl.get(args.process, ch) for ch in channels}
    if any(v is None for v in nominal.values()):
        raise SystemExit(f"nominal template of {args.process} missing in some channel")
    nom = np.concatenate([nominal[ch][0] for ch in channels])
    enom = np.sqrt(np.concatenate([nominal[ch][1] for ch in channels]))
    nbins = {ch: len(nominal[ch][0]) for ch in channels}
    print(f"{args.card}: {args.process}, channels {', '.join(f'{c} ({n} bins)' for c, n in nbins.items())}, "
          f"total yield {nom.sum():.1f}")
    floats = [p for p in card["params"] if len(p.split()) > 3 and p.split()[3] in (args.process, "*")]
    print("rateParams on it: " + ("; ".join(floats) if floats else "none"))

    rows = []
    for name, (typ, entries) in card_nuisances(card, args.process, args.groups).items():
        row = {"name": name, "type": typ}
        if typ == "lnN":
            vals = set(entries.values())
            down, up = lnn(next(iter(vals)))
            row.update({"rate_up": up, "rate_down": down, "rms_up": 0.0, "rms_down": 0.0,
                        "max_up": 0.0, "max_down": 0.0, "note": "" if len(vals) == 1 else "differs per channel"})
        else:
            per = {}
            for ud in ("Up", "Down"):
                hs = []
                for ch in channels:
                    h = tpl.get(args.process, ch, name + ud) if ch in entries else None
                    hs.append(h[0] if h is not None else nominal[ch][0])      # no entry: nominal
                per[ud] = shape_metrics(nom, enom, np.concatenate(hs))
            row.update({f"{k}_{ud.lower()}": per[ud][k] for ud in ("Up", "Down") for k in ("rate", "rms", "max", "chi2")})
            row.update({"ndf": per["Up"]["ndf"], "s_up": per["Up"]["s"], "s_down": per["Down"]["s"], "note": ""})
        rows.append(row)

    if args.hdamp:
        with open(args.hdamp) as f:
            hd = json.load(f)
        chs = list(hd["channels"])
        hn = np.concatenate([np.asarray(hd["channels"][c]["nominal"][0]) for c in chs])
        he = np.concatenate([np.asarray(hd["channels"][c]["nominal"][1]) for c in chs])
        row = {"name": "hdamp (4FS samples)", "type": "samples",
               "note": f"{hd['subsample']}, {'+'.join(hd['years'])}, sample stat. included"}
        for ud in ("up", "down"):
            hv = np.concatenate([np.asarray(hd["channels"][c][ud][0]) for c in chs])
            ev = np.concatenate([np.asarray(hd["channels"][c][ud][1]) for c in chs])
            m = shape_metrics(hn, he, hv, ev)
            row.update({f"rate_{ud}": m["rate"], f"rms_{ud}": m["rms"], f"rmssub_{ud}": m["rms_sub"],
                        f"max_{ud}": m["max"], f"chi2_{ud}": m["chi2"], "ndf": m["ndf"]})
            row["note"] += f"; {ud}: raw rms {100 * m['rms']:.1f}% vs {100 * m['rms_noise']:.1f}% expected from stat."
        rows.append(row)

    rows.sort(key=lambda r: -max(r.get("rmssub_up", r["rms_up"]), r.get("rmssub_down", r["rms_down"])))
    head = (f"{'nuisance':24} {'type':7} {'rate up':>8} {'rate dn':>8} {'shape rms up/dn [%]':>20} "
            f"{'max |s| up/dn [%]':>18} {'chi2 up/dn / ndf':>20}  note")
    print("\n" + head + "\n" + "-" * len(head))
    for r in rows:
        rms = (f"{100 * r.get('rmssub_up', r['rms_up']):.1f}/{100 * r.get('rmssub_down', r['rms_down']):.1f}"
               + (" *" if "rmssub_up" in r else ""))
        chi = (f"{r['chi2_up']:.0f}/{r['chi2_down']:.0f} / {r['ndf']}" if "chi2_up" in r else "-")
        mx = f"{100 * r['max_up']:.1f}/{100 * r['max_down']:.1f}"
        print(f"{r['name']:24} {r['type']:7} {r['rate_up']:8.3f} {r['rate_down']:8.3f} {rms:>20} {mx:>18} "
              f"{chi:>20}  {r['note']}")
    print("* stat-subtracted (independent samples); chi2 with the stat. of both samples. "
          "Rate parts are absorbed if the process floats.")

    with open(f"{args.outdir}/{args.process}_theory.csv", "w", newline="") as f:
        keys = ["name", "type", "rate_up", "rate_down", "rms_up", "rms_down", "rmssub_up", "rmssub_down",
                "max_up", "max_down", "chi2_up", "chi2_down", "ndf", "note"]
        wr = csv.DictWriter(f, keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)

    draw(rows, nbins, args)
    print(f"Wrote {args.outdir}/{args.process}_theory.csv, _theory_summary.png, _theory_bins.png")


def draw(rows, nbins, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    up_c, dn_c, hd_c = "#1f4e9c", "#d9730d", "#2a9d3c"
    names = [r["name"] for r in rows]
    y = np.arange(len(rows))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 0.38 * len(rows) + 1.6), sharey=True)
    for k, (ud, c) in enumerate((("up", up_c), ("down", dn_c))):
        off = (k - 0.5) * 0.36
        a1.barh(y + off, [100 * (r[f"rate_{ud}"] - 1) for r in rows], 0.34, color=c, label=ud)
        a2.barh(y + off, [100 * r.get(f"rmssub_{ud}", r[f"rms_{ud}"]) for r in rows], 0.34, color=c, label=ud)
    for i, r in enumerate(rows):
        if r["type"] == "samples":
            for a in (a1, a2):
                a.axhspan(i - 0.5, i + 0.5, color=hd_c, alpha=0.12, lw=0)
    a1.axvline(0, color="gray", lw=0.8)
    a1.set_xlabel("rate: variation / nominal - 1 [%]")
    a2.set_xlabel("shape only: yield-weighted RMS over the datacard bins [%]")
    a1.set_yticks(y, names, fontsize=8)
    a1.invert_yaxis()
    a1.legend(fontsize=8, frameon=False)
    for a in (a1, a2):
        a.grid(axis="x", alpha=0.3)
    fig.suptitle(f"{args.process}: theory nuisances in {os.path.basename(args.card)} "
                 f"(shaded: hdamp from the 4FS samples, stat-subtracted)", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{args.outdir}/{args.process}_theory_summary.{ext}")
    plt.close(fig)

    shapes = [r for r in rows if "s_up" in r][:args.top]
    if not shapes:
        return
    n = sum(nbins.values())
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    cmap = plt.get_cmap("tab10")
    for ax, ud in zip(axes, ("up", "down")):
        for k, r in enumerate(shapes):
            ax.stairs(100 * r[f"s_{ud}"], np.arange(n + 1), color=cmap(k), lw=1.4, label=r["name"])
        x0 = 0
        for ch, nb in nbins.items():
            ax.axvline(x0, color="gray", lw=0.6, ls=":")
            ax.text(x0 + 0.3, 0.97, ch, transform=ax.get_xaxis_transform(), va="top", fontsize=8)
            x0 += nb
        ax.axhline(0, color="gray", lw=0.8)
        ax.set_ylabel(f"{ud}: shape only [%]")
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=3, frameon=False, loc="lower left")
    axes[1].set_xlabel("datacard bin (channels concatenated)")
    fig.suptitle(f"{args.process}: shape-only effect of the largest card nuisances "
                 f"(variation scaled to the nominal total)", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{args.outdir}/{args.process}_theory_bins.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    main()
