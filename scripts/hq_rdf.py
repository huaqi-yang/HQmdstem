#!/usr/bin/env python3
"""HQmdstemkit: RDF plotting.

Usage:
  hq_rdf.py rdf FILE [--out out.png] [--pair 2]
  hq_rdf.py gr  FILE [--out out.png]
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, load_config, hq_home

USAGE = __doc__

def _load(f):
    with open(f) as fh:
        rows = []
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            rows.append([float(x) for x in s.split()])
    return np.array(rows)

def _plot(x, ys, labels, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.5, 3.2), dpi=200)
    for i, y in enumerate(ys):
        ax.plot(x, y, lw=1.5, label=labels[i] if labels else None)
    ax.set_xlabel("r ($\\AA$)")
    ax.set_ylabel("g(r)")
    ax.tick_params(direction="in", top=True, right=True)
    if labels:
        ax.legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"saved -> {out}")

def cmd_rdf(args):
    if len(args) < 1:
        fail(USAGE)
    f = args[0]
    out = "rdf.png"
    pair = None
    i = 1
    while i < len(args):
        if args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        elif args[i] == "--pair" and i+1 < len(args): pair = int(args[i+1]); i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    if not os.path.isfile(f):
        fail(f"file not found: {f}\n\n{USAGE}")
    data = _load(f)
    if data.shape[1] < 2:
        fail("rdf file needs at least 2 columns (r, g(r))\n\n" + USAGE)
    x = data[:, 0]
    if pair is not None:
        if pair >= data.shape[1]:
            fail(f"--pair {pair} out of range (columns={data.shape[1]})\n\n{USAGE}")
        ys = [data[:, pair]]
        labels = [f"pair {pair}"]
    else:
        ys = [data[:, j] for j in range(1, data.shape[1])]
        labels = [f"col {j}" for j in range(1, data.shape[1])]
    _plot(x, ys, labels, out)

def cmd_gr(args):
    if len(args) < 1:
        fail(USAGE)
    f = args[0]
    out = "gr.png"
    if "--out" in args:
        out = args[args.index("--out") + 1]
    if not os.path.isfile(f):
        fail(f"file not found: {f}\n\n{USAGE}")
    data = _load(f)
    x = data[:, 0]
    ys = [data[:, j] for j in range(1, data.shape[1])]
    _plot(x, ys, [f"col {j}" for j in range(1, data.shape[1])], out)


RDF41_COLORS = {"50K": "#4CAF50", "100K": "#8BC34A", "200K": "#FF9800", "250K": "#E53935"}
RDF41_TEMPS = ["50K", "100K", "200K", "250K"]
RDF41_STRUCTS = ["4", "8", "16", "32"]
RDF41_YLIMS = {"4": 9, "8": 13, "16": 13, "32": 13}

def _read_rdf_columns(filepath):
    r, g = [], []
    with open(filepath) as f:
        for line in f:
            s = line.split()
            if len(s) >= 3:
                try:
                    r.append(float(s[1]))
                    g.append(float(s[2]))
                except ValueError:
                    pass
    return np.array(r), np.array(g)

def cmd_4x1(args):
    """Port of RDF_4x1_merged.py: 4x1 RDF panels for SCRAPS 4/8/16/32."""
    base = load_config().get("rdf_base", os.path.join(hq_home(), "examples", "rdf_data"))
    out_dir = None
    i = 0
    while i < len(args):
        if args[i] == "--base" and i+1 < len(args):
            base = args[i+1]; i += 2
        elif args[i] == "--out" and i+1 < len(args):
            out_dir = args[i+1]; i += 2
        else:
            fail(f"unknown option {args[i]}\n\n{USAGE}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 19, "axes.linewidth": 2.6,
        "axes.labelsize": 19, "xtick.labelsize": 15, "ytick.labelsize": 15,
        "xtick.major.size": 7, "xtick.major.width": 2.6, "xtick.direction": "in",
        "ytick.major.size": 7, "ytick.major.width": 2.6, "ytick.direction": "in",
        "lines.linewidth": 3.6, "legend.fontsize": 15, "legend.frameon": False,
        "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })
    data = {}
    for s in RDF41_STRUCTS:
        data[s] = {}
        for temp in RDF41_TEMPS:
            path = os.path.join(base, s, temp, "1.txt")
            if os.path.isfile(path):
                data[s][temp] = _read_rdf_columns(path)
    fig, axes = plt.subplots(4, 1, figsize=(16/2.54, 32/2.54), sharex=True)
    fig.subplots_adjust(left=0.16, bottom=0.06, right=0.96, top=0.995, hspace=0.0)
    for i, (s, ax) in enumerate(zip(RDF41_STRUCTS, axes)):
        ax.grid(False)
        for temp in RDF41_TEMPS:
            if temp in data[s]:
                r, g = data[s][temp]
                g = g.copy(); g[r < 2.30] = 0.0
                ax.plot(r, g, color=RDF41_COLORS[temp], lw=2.3)
        ax.set_xlim(1.7, 3.2); ax.set_ylim(0, RDF41_YLIMS[s])
        ax.set_ylabel(r"$g(r)$", fontsize=17)
        if i < 3:
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel(r"$r$ ($\mathrm{\AA}$)")
        ax.tick_params(direction="in", top=True, right=True)
        for j, temp in enumerate(RDF41_TEMPS):
            ypos = 0.93 - j * 0.09
            ax.plot([0.78, 0.86], [ypos, ypos], transform=ax.transAxes,
                    color=RDF41_COLORS[temp], lw=2.8, solid_capstyle="butt", clip_on=False)
            ax.text(0.88, ypos, temp, transform=ax.transAxes,
                    color=RDF41_COLORS[temp], fontsize=14, ha="left", va="center")
    out_dir = out_dir or os.path.join(base, "RDF_figures")
    os.makedirs(out_dir, exist_ok=True)
    for fmt in ("png", "pdf", "svg"):
        out = os.path.join(out_dir, f"RDF_4x1.{fmt}")
        fig.savefig(out, format=fmt)
        print(f"[OK] {out}")
    plt.close(fig)

def main():
    if len(sys.argv) < 3:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "rdf": cmd_rdf(args)
    elif cmd == "gr": cmd_gr(args)
    elif cmd == "4x1": cmd_4x1(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()