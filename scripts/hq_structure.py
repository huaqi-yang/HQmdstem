#!/usr/bin/env python3
"""HQmdstemkit: ordered/disordered CuZn solid solution generator (fcc).

Usage:
  hq_structure.py ordered L12 N [--a 3.61] [--out L12.xyz]
  hq_structure.py ordered B2  N [--a 2.98] [--out B2.xyz]
  hq_structure.py ordered L10 N [--a 3.70] [--c 3.62] [--out L10.xyz]
  hq_structure.py disordered NCU NZN [--a 3.61] [--out random.xyz]
"""
import os, sys, random
import numpy as np

USAGE = __doc__

def write_xyz(out, lattice, species, pos):
    n = len(species)
    a, b, c = lattice
    with open(out, "w") as f:
        f.write(f"{n}\n")
        f.write(f'Lattice="{a[0]} {a[1]} {a[2]} {b[0]} {b[1]} {b[2]} {c[0]} {c[1]} {c[2]}" '
                f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for sp, p in zip(species, pos):
            f.write(f"{sp} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
    print(f"saved -> {out}")

def fcc_vectors(a):
    return (np.array([a, 0, 0]), np.array([0, a, 0]), np.array([0, 0, a]))

def cubic_sites(n):
    """n x n x n fcc supercell: conventional cell with 4 basis atoms."""
    basis = [(0,0,0),(0.5,0.5,0),(0.5,0,0.5),(0,0.5,0.5)]
    sites = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for b in basis:
                    sites.append((i+b[0], j+b[1], k+b[2]))
    return sites

def ordered_l12(n, a):
    """Cu3Zn L12: Zn on cube corners, Cu on face centers."""
    species, pos = [], []
    for s in cubic_sites(n):
        x, y, z = s
        corners = (abs(x-round(x)) < 1e-9 and abs(y-round(y)) < 1e-9 and abs(z-round(z)) < 1e-9)
        species.append("Zn" if corners else "Cu")
        pos.append(np.array([x, y, z]) * a)
    return fcc_vectors(a), species, np.array(pos)

def ordered_b2(n, a):
    """B2 CsCl-type: Cu at fcc sites, Zn at octahedral interstitial (bcc-like ordering)."""
    species, pos = [], []
    for s in cubic_sites(n):
        x, y, z = s
        site = (round((x % 1) * 2), round((y % 1) * 2), round((z % 1) * 2))
        species.append("Cu" if sum(site) % 2 == 0 else "Zn")
        pos.append(np.array([x, y, z]) * a)
    return fcc_vectors(a), species, np.array(pos)

def ordered_l10(n, a, c):
    """L10 CuZn: alternating Cu/Zn layers along z in fcc."""
    species, pos = [], []
    for s in cubic_sites(n):
        x, y, z = s
        layer = round(z) % 2
        species.append("Cu" if layer == 0 else "Zn")
        p = np.array([x, y, z * (c / a)])
        pos.append(p)
    a_vec = np.array([a, 0, 0]); b_vec = np.array([0, a, 0]); c_vec = np.array([0, 0, c])
    return (a_vec, b_vec, c_vec), species, np.array(pos)

def disordered(ncu, nzn, a):
    n = max(int(round((ncu + nzn) ** (1.0/3.0))), 2)
    sites = cubic_sites(n)
    random.shuffle(sites)
    total = ncu + nzn
    if len(sites) < total:
        n = int(np.ceil((total / 4.0) ** (1.0/3.0)))
        sites = cubic_sites(n)
    sel = sites[:total]
    species = ["Cu"] * ncu + ["Zn"] * nzn
    random.shuffle(species)
    pos = np.array(sel) * a
    return fcc_vectors(a), species, pos

def cmd_ordered(args):
    if len(args) < 2:
        print(USAGE); sys.exit(1)
    kind, n = args[0], int(args[1])
    a, c = 3.61, 3.62
    out = f"{kind}_{n}.xyz"
    i = 2
    while i < len(args):
        if args[i] == "--a" and i+1 < len(args): a = float(args[i+1]); i += 2
        elif args[i] == "--c" and i+1 < len(args): c = float(args[i+1]); i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: print(f"unknown option {args[i]}\n{USAGE}"); sys.exit(1)
    if kind == "L12": lat, sp, pos = ordered_l12(n, a)
    elif kind == "B2": lat, sp, pos = ordered_b2(n, a)
    elif kind == "L10": lat, sp, pos = ordered_l10(n, a, c)
    else: print(f"unknown ordered type {kind}\n{USAGE}"); sys.exit(1)
    write_xyz(out, lat, sp, pos)

def cmd_disordered(args):
    if len(args) < 2:
        print(USAGE); sys.exit(1)
    ncu, nzn = int(args[0]), int(args[1])
    a, out = 3.61, "CuZn_random.xyz"
    i = 2
    while i < len(args):
        if args[i] == "--a" and i+1 < len(args): a = float(args[i+1]); i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: print(f"unknown option {args[i]}\n{USAGE}"); sys.exit(1)
    lat, sp, pos = disordered(ncu, nzn, a)
    write_xyz(out, lat, sp, pos)

def main():
    if len(sys.argv) < 2:
        print(USAGE); sys.exit(1)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "ordered": cmd_ordered(args)
    elif cmd == "disordered": cmd_disordered(args)
    else: print(USAGE); sys.exit(1)

if __name__ == "__main__":
    main()