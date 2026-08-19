#!/usr/bin/env python3
"""HQmdstemkit: batch pretreatment (GPUMDkit 302/303 style).

Usage:
  hq_batch.py batch SRC [--prefix iter00] [--mode gpumd|lmp] [--out .]
"""
import os, sys, re, glob, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import read_extxyz_frames, load_config, hq_home, fail

USAGE = __doc__

def _copy(src, dst):
    if os.path.isfile(src):
        shutil.copyfile(src, dst)
        return True
    return False

def read_poscar(path):
    """Parse VASP POSCAR, tolerating missing element names, Cartesian/Direct,
    and Selective-dynamics tags (GPUMDkit-style vasp files)."""
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        raw = f.readlines()
    lines = [ln for ln in raw if ln.strip()]
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    elems = []
    counts = []
    i = 5
    while i < len(lines):
        parts = lines[i].split()
        try:
            counts = [int(x) for x in parts]
            break
        except ValueError:
            elems = parts
            i += 1
    if not counts:
        raise ValueError(f"cannot find atom counts in POSCAR: {path}")
    if not elems:
        elems = ["Cu", "Zn"][:len(counts)] if len(counts) <= 2 else [f"X{j+1}" for j in range(len(counts))]
    j = i + 1
    kind = "direct"
    if j < len(lines) and lines[j].strip().lower().startswith("selective"):
        j += 1
    if j < len(lines) and lines[j].strip().lower().startswith(("cartesian", "direct")):
        kind = "cartesian" if lines[j].strip().lower().startswith("c") else "direct"
        j += 1
    coords = []
    for ln in lines[j:j + sum(counts)]:
        coords.append([float(x) for x in ln.split()[:3]])
    if len(coords) < sum(counts):
        raise ValueError(f"POSCAR atom count mismatch: expected {sum(counts)}, got {len(coords)}")
    pos = np.array(coords)
    if kind == "direct":
        pos = pos @ lat
    species = []
    for el, n in zip(elems, counts):
        species += [el] * n
    return lat, species, pos

def write_model_xyz(path, lat, species, pos):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(species)}\n")
        f.write(f'Lattice="{lat[0,0]:.10f} {lat[0,1]:.10f} {lat[0,2]:.10f} '
                f'{lat[1,0]:.10f} {lat[1,1]:.10f} {lat[1,2]:.10f} '
                f'{lat[2,0]:.10f} {lat[2,1]:.10f} {lat[2,2]:.10f}" '
                f'Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for s, p in zip(species, pos):
            f.write(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")

def _copy_fn(mode, src):
    ex = hq_home()
    def copy_to(d):
        if mode == "gpumd":
            for f in ("run.in", "nep.txt"):
                if not _copy(os.path.join(src, f), os.path.join(d, f)):
                    _copy(os.path.join(ex, "examples", "gpumd", f), os.path.join(d, f))
        else:
            eam_dir = os.path.join(ex, "examples", "elastic", "eam")
            for f in ("run_wsl_lmp_10queue.py", "CuZn.eam.alloy", "lammps_elastic_eam.in",
                      "in.elastic", "init.mod", "potential.mod", "lammps.data", "presub.sh"):
                _copy(os.path.join(eam_dir, f), os.path.join(d, f))
    return copy_to

def cmd_batch(args):
    if len(args) < 1:
        fail(USAGE)
    src = args[0]
    prefix, mode, out = "iter00", "gpumd", "."
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--prefix" and i + 1 < len(args): prefix = args[i+1]; i += 2
        elif a == "--mode" and i + 1 < len(args): mode = args[i+1]; i += 2
        elif a == "--out" and i + 1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    os.makedirs(out, exist_ok=True)
    copy_fn = _copy_fn(mode, src)
    poscars = sorted(glob.glob(os.path.join(src, "POSCAR_*.vasp")))
    if not poscars:
        poscars = sorted(glob.glob(os.path.join(src, "*.vasp")))
    n_jobs = 0
    if poscars:
        for f in poscars:
            m = re.search(r"POSCAR_?(\d+)", os.path.basename(f))
            idx = m.group(1) if m else str(poscars.index(f) + 1)
            d = os.path.join(out, f"{prefix}_{idx}")
            os.makedirs(d, exist_ok=True)
            lat, species, pos = read_poscar(f)
            write_model_xyz(os.path.join(d, "model.xyz"), lat, species, pos)
            copy_fn(d)
            n_jobs += 1
            print(f"{os.path.basename(f)} -> {d}")
    else:
        xyz = os.path.join(src, "model.xyz")
        if not os.path.isfile(xyz):
            fail(f"no POSCAR_*.vasp or model.xyz in {src}\n\n{USAGE}")
        frames = read_extxyz_frames(xyz)
        for k, (n, header, atoms) in enumerate(frames, 1):
            d = os.path.join(out, f"{prefix}_{k}")
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "model.xyz"), "w", encoding="utf-8") as fh:
                fh.write(f"{n}\n{header}")
                for al in atoms:
                    fh.write(al if al.endswith("\n") else al + "\n")
            copy_fn(d)
            n_jobs += 1
            print(f"frame {k} -> {d}")
    print(f"batch done ({mode} mode): {n_jobs} jobs in {out}")

def main():
    if len(sys.argv) < 3:
        fail(USAGE)
    if sys.argv[1] == "batch":
        cmd_batch(sys.argv[2:])
    else:
        fail(USAGE)

if __name__ == "__main__":
    main()