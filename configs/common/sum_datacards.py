#!/usr/bin/env python3
"""
Sum per-year datacards (e.g. Run 2: 2016, 2017, 2018) into ONE datacard + shape file, for a single
fit over the integrated years instead of a combination of the per-year cards.

    cd configs/Run2/datacards
    python3 ../../common/sum_datacards.py datacard_nominal_disc_<model>_{2016,2017,2018}.txt \\
        -o datacard_nominal_disc_<model>_Run2

writes <output>.txt and <output>.root. The cards must have the same channels, processes and
binning (common SPANet binning). Per (channel, process):

  nominal                 sum of the years (variances add: autoMCStats sees the summed MC statistics)
  shape nuisances         Up = sum over the years of the year's Up, where a year without the
                          nuisance contributes its nominal. So a per-year nuisance (pu_2017) moves
                          only its own year; a correlated one (same name every year) moves all.
  lnN, same value in every year that has it and present in all years -> copied
  lnN in INTEGRATED_LNN   replaced by the integrated-years value (LUM recommendation)
  lnN in AVERAGED_LNN     average of the years' kappa_up and kappa_down (separately), weighted by
                          the process's total yield in each year (all channels)
  any other lnN that differs between years or misses in some -> error (no guessing)
  rateParam, autoMCStats  copied (union of the years' lines)
  groups                  union of the years' members
  data_obs                summed like the MC
"""
import argparse
import os
from collections import OrderedDict

import hist
import numpy as np
import uproot

# LUM POG, analyses using 2016+2017+2018 integrated
INTEGRATED_LNN = {
    "lumi_13TeV_1516": "1.0031",
    "lumi_13TeV_151617": "1.0018",
    "lumi_13TeV_15161718": "1.0064",
}
# lnN whose values legitimately differ slightly per year (computed per year from the same samples)
AVERAGED_LNN = ["ttx_qsc", "qqpdf"]


# ---------------------------------------------------------------- reading
def parse_card(path):
    """Datacard contents: shapes patterns, bins, columns [(bin, process, index)], nuisances
    {name: (type, [value per column])}, rateParam/autoMCStats lines, groups."""
    card = {"shapes": {}, "columns": None, "nuisances": OrderedDict(), "params": [], "autoMCStats": [],
            "groups": OrderedDict(), "dir": os.path.dirname(os.path.abspath(path))}
    lines = [l.split("#")[0].rstrip() for l in open(path)]
    lines = [l for l in lines if l.strip() and not set(l.strip()) <= {"-"}]
    start = next((i for i, l in enumerate(lines) if l.split()[0] == "imax"), None)
    if start is None:
        raise SystemExit(f"{path}: no imax line")
    lines = lines[start:]           # free-text header (e.g. "Datacard for 2016") is ignored
    for i, line in enumerate(lines):
        tok = line.split()
        key = tok[0]
        if key in ("imax", "jmax", "kmax", "observation"):
            continue
        if key == "shapes":
            card["shapes"][(tok[1], tok[2])] = tok[3:]
        elif key == "bin" and i + 1 < len(lines) and lines[i + 1].split()[0] == "process":
            names = lines[i + 1].split()[1:]
            indices = lines[i + 2].split()[1:]
            card["columns"] = list(zip(tok[1:], names, indices))
        elif key in ("bin", "process", "rate"):
            continue
        elif len(tok) > 2 and tok[1] == "group" and tok[2] == "=":
            card["groups"][tok[0]] = tok[3:]
        elif len(tok) > 1 and tok[1] == "rateParam":
            card["params"].append(" ".join(tok))
        elif len(tok) > 1 and tok[1] == "autoMCStats":
            card["autoMCStats"].append(" ".join(tok))
        elif card["columns"] and len(tok) == 2 + len(card["columns"]) and tok[1] in ("lnN", "shape"):
            card["nuisances"][tok[0]] = (tok[1], tok[2:])
        else:
            raise SystemExit(f"{path}: line not understood (not summed automatically): {line}")
    if card["columns"] is None:
        raise SystemExit(f"{path}: no bin/process block")
    return card


def shape_pattern(card, process, channel, systematic=False):
    """ROOT file and histogram name of a process (or its '<name>' systematic pattern) in a channel."""
    for key in ((process, channel), (process, "*"), ("*", channel), ("*", "*")):
        if key in card["shapes"]:
            spec = card["shapes"][key]
            break
    else:
        raise SystemExit(f"No shapes line for {process} in {channel}")
    name = spec[2] if systematic else spec[1]
    name = name.replace("$CHANNEL", channel).replace("$PROCESS", process)
    return os.path.join(card["dir"], spec[0]), name


class Templates:
    """Histogram access for one card (files opened once)."""
    def __init__(self, card):
        self.card, self.files = card, {}

    def get(self, process, channel, systematic=None):
        path, name = shape_pattern(self.card, process, channel, systematic is not None)
        if systematic is not None:
            name = name.replace("$SYSTEMATIC", systematic)
        f = self.files.setdefault(path, uproot.open(path))
        if name not in f:
            return None
        h = f[name]
        return h.values(), h.variances(), h.axis().edges()


def lnn(value):
    """(kappa_down, kappa_up) of an lnN entry ('1.05' or '0.95/1.06')."""
    if "/" in value:
        down, up = (float(v) for v in value.split("/"))
        return down, up
    k = float(value)
    return 1 / k, k


# ---------------------------------------------------------------- summing
def check_same_structure(cards, paths):
    ref = [(b, p, i) for b, p, i in cards[0]["columns"]]
    for card, path in zip(cards[1:], paths[1:]):
        if card["columns"] != ref:
            raise SystemExit(f"{path}: bins / processes differ from {paths[0]}")


def process_yields(cards, templates):
    """{(year index, process): total nominal yield over all channels}."""
    out = {}
    for y, (card, tpl) in enumerate(zip(cards, templates)):
        for b, p, _ in card["columns"]:
            values = tpl.get(p, b)[0]
            out[(y, p)] = out.get((y, p), 0.0) + float(np.sum(values))
    return out


def summed_lnn(name, cards, paths, yields):
    """Value per column of an lnN nuisance in the summed card."""
    columns = cards[0]["columns"]
    per_year = [card["nuisances"].get(name) for card in cards]
    present = [y for y, n in enumerate(per_year) if n is not None]
    out = []
    for c, (b, p, _) in enumerate(columns):
        values = {y: per_year[y][1][c] for y in present if per_year[y][1][c] != "-"}
        if not values:
            out.append("-")
        elif name in INTEGRATED_LNN:
            out.append(INTEGRATED_LNN[name])
        elif name in AVERAGED_LNN:
            weights = np.array([yields[(y, p)] for y in values])
            kd = np.average([lnn(v)[0] for v in values.values()], weights=weights)
            ku = np.average([lnn(v)[1] for v in values.values()], weights=weights)
            out.append(f"{kd:.4f}/{ku:.4f}")
        else:
            if len(present) != len(cards) or len(set(values.values())) != 1 or len(values) != len(present):
                detail = {os.path.basename(paths[y]): per_year[y][1][c] if per_year[y] else "absent"
                          for y in range(len(cards))}
                raise SystemExit(f"lnN '{name}' differs between the years for {b}/{p}: {detail}\n"
                                 f"Add it to INTEGRATED_LNN or AVERAGED_LNN in sum_datacards.py if intended.")
            out.append(next(iter(values.values())))
    return out


def summed_shape(name, cards, templates, channel, process, column):
    """(Up, Down) summed over the years; years without the nuisance contribute their nominal."""
    up = down = None
    used = False
    for card, tpl in zip(cards, templates):
        nom = tpl.get(process, channel)[0]
        n = card["nuisances"].get(name)
        if n is not None and n[1][column] != "-":
            if n[1][column] != "1":
                raise SystemExit(f"shape '{name}' has scale {n[1][column]} for {channel}/{process}: only 1 supported")
            u, d = tpl.get(process, channel, f"{name}Up"), tpl.get(process, channel, f"{name}Down")
            if u is None or d is None:
                raise SystemExit(f"Missing {name}Up/Down template for {channel}/{process} in {card['dir']}")
            u, d, used = u[0], d[0], True
        else:
            u = d = nom
        up = u if up is None else up + u
        down = d if down is None else down + d
    return (up, down) if used else None


def to_hist(values, variances, edges):
    h = hist.Hist(hist.axis.Variable(edges), storage=hist.storage.Weight())
    h.view().value = values
    h.view().variance = variances
    return h


# ---------------------------------------------------------------- writing
def write(cards, paths, out):
    templates = [Templates(c) for c in cards]
    columns = cards[0]["columns"]
    channels = list(OrderedDict.fromkeys(b for b, _, _ in columns))
    yields = process_yields(cards, templates)
    root_name = os.path.basename(out) + ".root"
    hists = {}

    def add(name, values, variances, edges):
        hists[name] = to_hist(values, variances, edges)

    # nominal and data
    for c, (b, p, _) in enumerate(columns):
        parts = [tpl.get(p, b) for tpl in templates]
        add(f"{b}_{p}", sum(x[0] for x in parts), sum(x[1] for x in parts), parts[0][2])
    for b in channels:
        parts = [tpl.get("data_obs", b) for tpl in templates]
        if all(x is not None for x in parts):
            add(f"{b}_data_obs", sum(x[0] for x in parts), sum(x[1] for x in parts), parts[0][2])

    # nuisances (union, in order of first appearance)
    names = list(OrderedDict.fromkeys(n for card in cards for n in card["nuisances"]))
    nuisance_lines = []
    for name in names:
        kinds = {card["nuisances"][name][0] for card in cards if name in card["nuisances"]}
        if len(kinds) != 1:
            raise SystemExit(f"Nuisance '{name}' has different types in the cards: {kinds}")
        kind = kinds.pop()
        if kind == "lnN":
            values = summed_lnn(name, cards, paths, yields)
        else:
            values = []
            for c, (b, p, _) in enumerate(columns):
                summed = summed_shape(name, cards, templates, b, p, c)
                if summed is None:
                    values.append("-")
                    continue
                nom_edges = templates[0].get(p, b)[2]
                nom_var = hists[f"{b}_{p}"].view().variance
                add(f"{b}_{p}_{name}Up", summed[0], nom_var, nom_edges)
                add(f"{b}_{p}_{name}Down", summed[1], nom_var, nom_edges)
                values.append("1")
        if any(v != "-" for v in values):
            nuisance_lines.append((name, kind, values))

    with uproot.recreate(out + ".root") as f:
        for name, h in hists.items():
            f[name] = h

    w = max(len(x) for x in [b for b, _, _ in columns] + [p for _, p, _ in columns]
            + [v for _, _, vals in nuisance_lines for v in vals]) + 2
    head = max(len(x) for x in [n for n, _, _ in nuisance_lines] + ["observation"]) + 8
    col = lambda vals: "".join(f"{v:<{w}}" for v in vals)
    lines = [f"# Sum of: {', '.join(os.path.basename(p) for p in paths)} (sum_datacards.py)",
             "imax * number of bins", "jmax * number of processes minus 1",
             "kmax * number of nuisance parameters", "-" * 80,
             f"shapes data_obs * {root_name} $CHANNEL_data_obs",
             f"shapes * * {root_name} $CHANNEL_$PROCESS $CHANNEL_$PROCESS_$SYSTEMATIC", "-" * 80,
             f"{'bin':<{head}}" + col(channels), f"{'observation':<{head}}" + col(["-1"] * len(channels)), "-" * 80,
             f"{'bin':<{head}}" + col([b for b, _, _ in columns]),
             f"{'process':<{head}}" + col([p for _, p, _ in columns]),
             f"{'process':<{head}}" + col([i for _, _, i in columns]),
             f"{'rate':<{head}}" + col(["-1"] * len(columns)), "-" * 80]
    for name, kind, values in nuisance_lines:
        lines.append(f"{name:<{head - 7}} {kind:<6}" + col(values))
    lines += list(OrderedDict.fromkeys(l for card in cards for l in card["params"]))
    lines += list(OrderedDict.fromkeys(l for card in cards for l in card["autoMCStats"]))
    groups = OrderedDict()
    for card in cards:
        for g, members in card["groups"].items():
            groups.setdefault(g, [])
            groups[g] += [m for m in members if m not in groups[g]]
    lines += [f"{g} group = {' '.join(m)}" for g, m in groups.items()]
    with open(out + ".txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote {out}.txt and {out}.root: {len(channels)} channels, {len(columns)} columns, "
          f"{len(nuisance_lines)} nuisances")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cards", nargs="+", help="per-year datacard .txt files")
    ap.add_argument("-o", "--output", required=True, help="output path without extension")
    args = ap.parse_args()
    cards = [parse_card(p) for p in args.cards]
    check_same_structure(cards, args.cards)
    write(cards, args.cards, args.output)


if __name__ == "__main__":
    main()
