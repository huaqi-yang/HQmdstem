#!/usr/bin/env python3
"""HQmdstemkit: microstructure auto-analysis.

Commands:
  hq_micro.py chain IMAGE [options]        # white-chain length statistics (TEM)
  hq_micro.py cluster XYZ [--rcut R] [--out cluster.csv]
  hq_micro.py segregate XYZ [--axis z] [--bins 20] [--out seg.csv]
  hq_micro.py twin XYZ [--axis z] [--tol 0.3]
  hq_micro.py grain XYZ [--rcut R] [--out grain.csv]
  hq_micro.py orientation XYZ [--out orient.png]

Notes:
  cluster/grain use a simple neighbor-graph (no PBC) approximation;
  orientation uses nearest-neighbour bond angles; twin uses layer-stacking
  composition analysis. For production use OVITO/CNA or a dedicated tool.
"""
import os, sys, math
import numpy as np
from collections import deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import read_extxyz_frames, fail

USAGE = __doc__
CHAIN_CORE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hq_tem_core.py")
CHAIN_DEFAULT = CHAIN_CORE

def load_frame(path, frame_index=0):
    frames = read_extxyz_frames(path)
    if not frames:
        fail(f"no frames in {path}\n\n{USAGE}")
    n, header, atoms = frames[frame_index]
    sp, pos = [], []
    for al in atoms:
        c = al.split()
        if len(c) >= 4:
            sp.append(c[0])
            pos.append([float(c[1]), float(c[2]), float(c[3])])
    return sp, np.array(pos), header

def estimate_rcut(pos):
    d = []
    for i in range(min(len(pos), 2000)):
        dist = np.linalg.norm(pos - pos[i], axis=1)
        dist[dist < 1e-6] = np.inf
        d.append(dist.min())
    return float(np.median(d)) * 1.3

def neighbor_components(pos, rcut):
    n = len(pos)
    adj = [[] for _ in range(n)]
    for i in range(n):
        dv = pos - pos[i]
        dist2 = np.einsum("ij,ij->i", dv, dv)
        for j in np.where(dist2 < rcut**2)[0]:
            if j != i:
                adj[i].append(int(j))
    seen = np.zeros(n, bool)
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        q = deque([i]); seen[i] = True; comp = []
        while q:
            u = q.popleft(); comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True; q.append(v)
        comps.append(comp)
    return comps

def cmd_chain(args):
    if len(args) < 1:
        fail(USAGE)
    import subprocess
    core = CHAIN_CORE if os.path.isfile(CHAIN_CORE) else CHAIN_DEFAULT
    if not os.path.isfile(core):
        fail(f"chain analysis core not found\n\n{USAGE}")
    subprocess.run([sys.executable, core, args[0]] + args[1:], check=False)

def cmd_cluster(args):
    if len(args) < 1:
        fail(USAGE)
    xyz = args[0]
    rcut = None; out = "cluster.csv"
    i = 1
    while i < len(args):
        if args[i] == "--rcut" and i+1 < len(args): rcut = float(args[i+1]); i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    sp, pos, _ = load_frame(xyz)
    if rcut is None:
        rcut = estimate_rcut(pos)
    comps = neighbor_components(pos, rcut)
    sizes = np.array([len(c) for c in comps])
    sizes.sort()
    print(f"atoms={len(pos)}  rcut={rcut:.3f} A  clusters={len(comps)}")
    print(f"cluster size: min={sizes.min()} median={np.median(sizes):.0f} max={sizes.max()}")
    with open(out, "w") as f:
        f.write("cluster_id,size\n")
        for k, s in enumerate(sizes, 1):
            f.write(f"{k},{s}\n")
    print(f"saved -> {out}")

def cmd_segregate(args):
    if len(args) < 1:
        fail(USAGE)
    xyz = args[0]
    axis = "z"; bins = 20; out = "segregation.csv"
    i = 1
    while i < len(args):
        if args[i] == "--axis" and i+1 < len(args): axis = args[i+1]; i += 2
        elif args[i] == "--bins" and i+1 < len(args): bins = int(args[i+1]); i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    sp, pos, _ = load_frame(xyz)
    idx = {"x": 0, "y": 1, "z": 2}[axis.lower()]
    c = pos[:, idx]
    edges = np.linspace(c.min(), c.max(), bins + 1)
    with open(out, "w") as f:
        f.write(f"bin_center,{axis},n_total,n_Cu,n_Zn,x_Zn\n")
        for k in range(bins):
            m = (c >= edges[k]) & (c < edges[k+1])
            n_cu = sum(1 for s, mm in zip(sp, m) if mm and s == "Cu")
            n_zn = sum(1 for s, mm in zip(sp, m) if mm and s == "Zn")
            n = n_cu + n_zn
            x_zn = n_zn / n if n else 0.0
            f.write(f"{(edges[k]+edges[k+1])/2:.4f},{axis},{n},{n_cu},{n_zn},{x_zn:.4f}\n")
    print(f"saved -> {out}")

def cmd_twin(args):
    if len(args) < 1:
        fail(USAGE)
    xyz = args[0]
    axis = "z"; tol = 0.3
    i = 1
    while i < len(args):
        if args[i] == "--axis" and i+1 < len(args): axis = args[i+1]; i += 2
        elif args[i] == "--tol" and i+1 < len(args): tol = float(args[i+1]); i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    sp, pos, _ = load_frame(xyz)
    idx = {"x": 0, "y": 1, "z": 2}[axis.lower()]
    order = np.argsort(pos[:, idx])
    c = pos[order, idx]
    sp_ord = [sp[k] for k in order]
    layers = []
    cur = [order[0]]
    for k in range(1, len(order)):
        if c[k] - c[k-1] > tol:
            layers.append(cur); cur = []
        cur.append(order[k])
    layers.append(cur)
    print(f"layers along {axis}: {len(layers)}")
    # twin heuristic: layer composition alternation reverses at a twin boundary
    xzn = []
    for lay in layers:
        nz = sum(1 for k in lay if sp[k] == "Zn")
        xzn.append(nz / len(lay))
    xzn = np.array(xzn)
    flip = np.sign(np.diff(xzn))
    flips = np.where(np.diff(flip) != 0)[0]
    if len(flips):
        print(f"twin-boundary candidates at layer {[int(f+1) for f in flips]} "
              f"(composition pattern reversal)")
    else:
        print("no clear twin-boundary candidate detected")

def cmd_grain(args):
    if len(args) < 1:
        fail(USAGE)
    xyz = args[0]
    rcut = None; out = "grain.csv"
    i = 1
    while i < len(args):
        if args[i] == "--rcut" and i+1 < len(args): rcut = float(args[i+1]); i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    sp, pos, _ = load_frame(xyz)
    if rcut is None:
        rcut = estimate_rcut(pos)
    comps = neighbor_components(pos, rcut)
    sizes = np.array([len(c) for c in comps])
    sizes.sort()
    print(f"grain-like domains (neighbor-graph approx): {len(comps)}, "
          f"median size {np.median(sizes):.0f} atoms")
    with open(out, "w") as f:
        f.write("domain_id,size_atoms\n")
        for k, s in enumerate(sizes, 1):
            f.write(f"{k},{s}\n")
    print(f"saved -> {out}")

def cmd_orientation(args):
    if len(args) < 1:
        fail(USAGE)
    xyz = args[0]
    out = "orientation.png"
    if "--out" in args:
        out = args[args.index("--out") + 1]
    sp, pos, _ = load_frame(xyz)
    angles = []
    for i in range(min(len(pos), 4000)):
        dv = pos - pos[i]
        dist2 = np.einsum("ij,ij->i", dv, dv)
        dist2[dist2 < 1e-8] = np.inf
        nb = np.argsort(dist2)[:2]
        if len(nb) >= 2:
            a, b = dv[nb[0]], dv[nb[1]]
            cosv = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
            angles.append(math.degrees(math.acos(max(-1, min(1, cosv)))))
    angles = np.array(angles)
    print(f"local NN bond angle: mean={angles.mean():.2f} deg, "
          f"median={np.median(angles):.2f} deg, std={angles.std():.2f} deg")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.2), dpi=200)
        ax.hist(angles, bins=36, range=(0, 180), color="#4CAF50", edgecolor="black", alpha=0.8)
        ax.set_xlabel("Nearest-neighbour bond angle (deg)")
        ax.set_ylabel("Count")
        ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout(); fig.savefig(out, dpi=200)
        print(f"saved -> {out}")
    except Exception as e:
        print(f"[warn] plotting skipped: {e}")

def main():
    if len(sys.argv) < 3:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "chain": cmd_chain(args)
    elif cmd == "cluster": cmd_cluster(args)
    elif cmd == "segregate": cmd_segregate(args)
    elif cmd == "twin": cmd_twin(args)
    elif cmd == "grain": cmd_grain(args)
    elif cmd == "orientation": cmd_orientation(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()