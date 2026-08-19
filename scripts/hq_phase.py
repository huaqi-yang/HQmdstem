#!/usr/bin/env python3
"""HQmdstemkit: Cu-Zn convex hull and experimental phase diagram.

Usage:
  hq_phase.py hull [--csv phases.csv] [--out convex_hull.png]
  hq_phase.py exp  [--csv exp_CuZn_boundaries.csv] [--out exp_phase.png]

CSV example (hull): phase,x_Zn,e_rel(eV/atom)
  Cu,0,0
  Cu3Zn-L12,0.25,-0.07
  CuZn-B2,0.5,-0.12
  Zn,1,0
"""
import os, sys, csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail

USAGE = __doc__

DEFAULT_HULL_DATA = [
    ("alpha-Cu FCC", 0.00, 0.0000), ("alpha-Cu3Zn L12", 0.25, -0.0366),
    ("beta-CuZn B2", 0.50, -0.0806), ("beta'-CuZn B2", 0.50, -0.0792),
    ("gamma-Cu5Zn8", 0.615, -0.0267), ("delta-CuZn3", 0.75, -0.0115),
    ("epsilon-CuZn4", 0.813, -0.0070), ("eta-Zn HCP", 1.00, 0.0000),
]

def monotone_chain(points):
    pts = sorted(points)
    if len(pts) <= 2:
        return pts
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def cmd_hull(args):
    data = DEFAULT_HULL_DATA
    csv_path = None
    out = "CuZn_convex_hull.png"
    i = 0
    while i < len(args):
        if args[i] == "--csv" and i+1 < len(args): csv_path = args[i+1]; i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    if csv_path:
        if not os.path.isfile(csv_path):
            fail(f"CSV not found: {csv_path}\n\n{USAGE}")
        data = []
        with open(csv_path, encoding="utf-8-sig", errors="replace") as f:
            for row in csv.DictReader(f):
                data.append((row["phase"], float(row["x_zn"]), float(row["e_rel"])))
    pts = [(x, e) for _, x, e in data]
    hull = monotone_chain(pts)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5), dpi=200)
    hx = [p[0] for p in hull]; hy = [p[1] for p in hull]
    ax.plot(hx, hy, "-", color="#0F4D92", lw=2)
    for name, x, e in data:
        ax.scatter(x, e, s=45, facecolors="none", edgecolors="#B64342", linewidths=1.5)
        if x < 0.95:
            ax.annotate(name, (x, e), textcoords="offset points", xytext=(5, 5), fontsize=7)
    ax.set_xlabel("Zn fraction $x_{Zn}$")
    ax.set_ylabel("Relative energy (eV/atom)")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"saved -> {out}")

DEFAULT_BOUNDARIES = {
    "alpha_ab": [(25,30.0),(200,31.0),(400,32.5),(454,36.8),(500,37.2),(558,37.3),(600,36.8),(700,35.8),(800,34.5),(902,32.1)],
    "beta_left": [(454,44.8),(500,44.5),(558,44.2),(600,43.8),(700,43.2),(800,42.5),(834,41.5)],
    "beta_right": [(454,48.5),(500,51.0),(558,55.0),(600,56.0),(700,57.0),(800,57.5),(834,57.8)],
    "gamma_left": [(454,48.5),(500,55.0),(558,57.6),(600,58.5),(700,58.8),(800,60.0),(834,60.5)],
    "gamma_right": [(558,68.0),(600,67.5),(700,66.5),(800,65.0)],
    "epsilon_left": [(558,73.0),(600,74.0),(700,75.5),(800,77.0)],
    "epsilon_right": [(419,85.0),(500,83.5),(558,83.0),(600,82.5),(700,81.5)],
    "eta_left": [(419,98.3),(500,98.5),(558,98.5),(600,99.0),(700,99.0)],
}

def cmd_exp(args):
    boundaries = DEFAULT_BOUNDARIES
    csv_path = None
    out = "exp_CuZn_phase_diagram.png"
    i = 0
    while i < len(args):
        if args[i] == "--csv" and i+1 < len(args): csv_path = args[i+1]; i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    if csv_path and os.path.isfile(csv_path):
        boundaries = {}
        with open(csv_path, encoding="utf-8-sig", errors="replace") as f:
            for row in csv.DictReader(f):
                boundaries.setdefault(row["boundary"], []).append(
                    (float(row["T_C"]), float(row["x_Zn_at_pct"])))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=200)
    for b, pts in boundaries.items():
        p = np.array(pts)
        ax.plot(p[:, 1], p[:, 0], lw=1.8, label=b)
    ax.set_xlim(0, 100); ax.set_ylim(0, 1000)
    ax.set_xlabel("Zn (at.%)"); ax.set_ylabel("Temperature ($^\\circ$C)")
    ax.set_title("Cu-Zn Phase Diagram (Experimental)")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(True, alpha=0.15)
    fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"saved -> {out}")


# NEP89 / NEPCuZn / DFT (ABACUS) ordered + SRO energies from the paper figure script
HULL_DFT_DATA = {
    "NEP89": {
        "ref": (-4.2760, -1.4862),
        "ordered": [("Cu3Zn L12", 3, 1, -3.6498), ("CuZn B2", 1, 1, -2.9709),
                    ("CuZn B2p", 1, 1, -2.9762), ("Cu5Zn8", 20, 32, -2.6312),
                    ("CuZn3", 13, 39, -2.2529), ("CuZn4", 3, 13, -2.0629)],
        "sro": [("SRO-8", 4, 4, -2.9691), ("SRO-16", 8, 8, -2.9701), ("SRO-32", 16, 16, -2.9701)],
    },
    "NEPCuZn": {
        "ref": (-3.5382, -1.2394),
        "ordered": [("Cu3Zn L12", 3, 1, -3.0385), ("CuZn B2", 1, 1, -2.5022),
                    ("CuZn B2p", 1, 1, -2.5010), ("Cu5Zn8", 20, 32, -2.2307),
                    ("CuZn3", 13, 39, -1.8725), ("CuZn4", 3, 13, -1.7299)],
        "sro": [("SRO-8", 4, 4, -2.4047), ("SRO-16", 8, 8, -2.3952), ("SRO-32", 16, 16, -2.3960)],
    },
    "DFT": {
        "ref": (-3.5506, -1.2690),
        "ordered": [("Cu3Zn L12", 3, 1, -3.0440), ("CuZn B2", 1, 1, -2.4964),
                    ("CuZn B2p", 1, 1, -2.4949), ("CuZn4", 3, 13, -1.7501)],
        "sro": [("SRO-4", 1, 1, -2.4169), ("SRO-8", 4, 4, -2.4016),
                ("SRO-16", 8, 8, -2.3927), ("SRO-32", 16, 16, -2.3937)],
    },
}

def cmd_hull_dft(args):
    """Port of convex_hull_with_dft.py: NEP89 + NEPCuZn + DFT convex hull."""
    out = "CuZn_convex_hull_with_dft.png"
    i = 0
    while i < len(args):
        if args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    def form_energy(e, n_cu, n_zn, refs):
        e_cu, e_zn = refs
        n = n_cu + n_zn
        return e - (n_cu/n)*e_cu - (n_zn/n)*e_zn
    model_colors = {"NEP89": "#0F4D92", "NEPCuZn": "#B64342", "DFT": "#2ECC71"}
    model_markers = {"NEP89": "o", "NEPCuZn": "D", "DFT": "^"}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=200)
    for model, cfg in HULL_DFT_DATA.items():
        pts = [(0.0, 0.0), (1.0, 0.0)]
        for name, n_cu, n_zn, e in cfg["ordered"]:
            x = n_zn/(n_cu+n_zn)
            df = form_energy(e, n_cu, n_zn, cfg["ref"])
            if abs(df) < 0.5:
                pts.append((x, df))
        hull = monotone_chain(pts)
        hx = [p[0] for p in hull]; hy = [p[1] for p in hull]
        ax.plot(hx, hy, "--", color=model_colors[model], lw=1.8, alpha=0.8)
        for name, n_cu, n_zn, e in cfg["ordered"]:
            x = n_zn/(n_cu+n_zn)
            df = form_energy(e, n_cu, n_zn, cfg["ref"])
            if abs(df) < 0.5:
                ax.scatter(x, df, s=55, facecolors="none", edgecolors=model_colors[model],
                           linewidths=1.4, marker=model_markers[model])
        for name, n_cu, n_zn, e in cfg["sro"]:
            x = n_zn/(n_cu+n_zn)
            df = form_energy(e, n_cu, n_zn, cfg["ref"])
            ax.scatter(x, df, s=45, facecolors="none", edgecolors=model_colors[model],
                       linewidths=1.6, marker="s", alpha=0.8)
    ax.axhline(0, color="#272727", lw=0.6, alpha=0.3)
    ax.set_xlabel("Zn fraction $x_{Zn}$")
    ax.set_ylabel("Relative total energy (eV/atom)")
    ax.set_xlim(-0.03, 1.03)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"saved -> {out}")

def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "hull": cmd_hull(args)
    elif cmd == "hull-dft": cmd_hull_dft(args)
    elif cmd == "exp": cmd_exp(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()