#!/usr/bin/env python3
"""HQmdstemkit: ABACUS SCF batch pretreatment / extract / energy shift.

Usage:
  hq_abacus.py pretreat SRC [--prefix iter00] [--out OUT] [--pp DIR] [--orb DIR] [--func LDA|PBE|PBE-D3]
  hq_abacus.py extract DIR [--out total_train.xyz]
  hq_abacus.py shift IN.xyz OUT.xyz
  hq_abacus.py sbatch [DIR]
"""
import os, sys, glob, re, shutil, subprocess
import numpy as np

try:
    from ase.io import read as _ase_read
    from ase.io import write as _ase_write
    HAVE_ASE = True
except Exception:
    HAVE_ASE = False

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_batch import read_poscar
from hq_common import fail, hq_home, read_extxyz_frames

USAGE = __doc__
EX_DIR = os.path.join(hq_home(), "examples", "abacus")
PP_DIR = os.path.join(EX_DIR, "apns-pseudopotentials-v1")
ORB_DIR = os.path.join(EX_DIR, "apns-orbitals-precision-v1")

ABACUS_INPUT = """INPUT_PARAMETERS
pseudo_dir      {pseudo}
orbital_dir     {orbital}
calculation     scf
basis_type      lcao
ecutwfc         100
scf_thr         1e-7
scf_nmax        100
device          gpu
ks_solver       cusolver
precision       double
gamma_only      0
kspacing        0.14
smearing_method gauss
smearing_sigma  0.02
mixing_type     pulay
mixing_beta     0.6
cal_force       1
cal_stress      1
out_charge      1
dft_functional  {func}
{vdw}"""


ABACUS_SBATCH = """#!/bin/bash
#SBATCH --job-name=CuZn-scf
#SBATCH --partition=4V100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --qos=flood-1o2gpu
#SBATCH --array=1-{n}%50

module load abacus/LTSv3.10.1-sm70-auto
WORK_DIR="$PWD/iter00_${{SLURM_ARRAY_TASK_ID}}"
cd "$WORK_DIR"

abacus > scf.output 2>&1
"""
def _stru_header(species):
    """Header used by the battle-tested SAI prepare_all2-8.sh (ABACUS STRU)."""
    cnt = {}
    for s in species:
        cnt[s] = cnt.get(s, 0) + 1
    masses = {"Cu": 63.546, "Zn": 65.38}
    lines = ["ATOMIC_SPECIES"]
    for el in sorted(cnt):
        pp_cands = sorted(glob.glob(os.path.join(PP_DIR, el + "_*.upf")) +
                          glob.glob(os.path.join(PP_DIR, el + "*.UPF")) +
                          glob.glob(os.path.join(PP_DIR, el + ".upf")))
        pp = os.path.basename(pp_cands[0]) if pp_cands else f"{el}.upf"
        lines.append(f"{el} {masses.get(el, 1.0)} {pp}")
    lines.append("NUMERICAL_ORBITAL")
    for el in sorted(cnt):
        orb_cands = sorted(glob.glob(os.path.join(ORB_DIR, el + "_*.orb")))
        orb = os.path.basename(orb_cands[0]) if orb_cands else f"{el}.orb"
        lines.append(orb)
    return "\n".join(lines) + "\n"


def write_stru(path, lat, species, pos):
    """Write STRU in the same format as abacustest/prepare_all2-8.sh output."""
    cnt = {}
    for s in species:
        cnt[s] = cnt.get(s, 0) + 1
    with open(path, "w", encoding="utf-8") as f:
        f.write(_stru_header(species))
        f.write("LATTICE_CONSTANT\n1.889726\n")
        f.write("LATTICE_VECTORS\n")
        for row in lat:
            f.write(f"    {row[0]:.10f} {row[1]:.10f} {row[2]:.10f}\n")
        f.write("ATOMIC_POSITIONS\nCartesian\n")
        for el in sorted(cnt):
            f.write(f"{el}\n0.000000\n{cnt[el]}\n")
            for s, p in zip(species, pos):
                if s == el:
                    f.write(f"    {p[0]:.10f} {p[1]:.10f} {p[2]:.10f} 1 1 1\n")


def _write_poscar_vasp(path, lat, species, pos):
    cnt = {}
    for s in species:
        cnt[s] = cnt.get(s, 0) + 1
    elems = sorted(cnt)
    with open(path, "w", encoding="utf-8") as f:
        f.write("converted by HQmdstemkit\n1.0\n")
        for row in lat:
            f.write(f"    {row[0]:.12f} {row[1]:.12f} {row[2]:.12f}\n")
        f.write(" ".join(elems) + "\n")
        f.write(" ".join(str(cnt[e]) for e in elems) + "\n")
        f.write("Cartesian\n")
        for e, p in zip(species, pos):
            f.write(f"{e} {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}\n")



def _cell_from_header(header):
    m = re.search(r"Lattice\s*=\s*[\"']?([-+0-9.eE\s]+)[\"']?", header)
    if m:
        vals = m.group(1).split()
        if len(vals) >= 9:
            return np.array([[float(x) for x in vals[i*3:i*3+3]] for i in range(3)])
    return None


def _frames_to_species_pos(fr):
    n = fr[0]
    sp, pos = [], []
    for al in fr[2:2 + n]:
        c = al.split()
        if len(c) >= 4:
            sp.append(c[0].capitalize())
            pos.append([float(c[1]), float(c[2]), float(c[3])])
    return sp, np.array(pos)
def _convert_to_poscars(src, out):
    """Auto-detect POSCAR_*.vasp / POSCAR / *.xyz / *.data and make POSCAR_1..n.vasp."""
    os.makedirs(out, exist_ok=True)
    existing = sorted(glob.glob(os.path.join(src, "POSCAR_*.vasp")))
    next_idx = 1
    for f in existing:
        m = re.search(r"POSCAR_?(\d+)", os.path.basename(f))
        idx = int(m.group(1)) if m else next_idx
        next_idx = max(next_idx, idx + 1)
        dst = os.path.join(out, os.path.basename(f))
        if os.path.abspath(f) != os.path.abspath(dst):
            shutil.copyfile(f, dst)
    single = os.path.join(src, "POSCAR")
    if os.path.isfile(single):
        dst = os.path.join(out, f"POSCAR_{next_idx}.vasp")
        if os.path.abspath(single) != os.path.abspath(dst):
            shutil.copyfile(single, dst)
        print(f"[ABACUS] POSCAR -> {os.path.basename(dst)}")
        next_idx += 1
    for f in sorted(glob.glob(os.path.join(src, "*.xyz"))):
        base = os.path.basename(f)
        conv = []
        if HAVE_ASE:
            try:
                atoms_list = _ase_read(f, index=":")
                conv = [("ase", a) for a in atoms_list]
            except Exception as e:
                print(f"[WARN] ase read failed {f}: {e}")
        if not conv:
            try:
                for fr in read_extxyz_frames(f):
                    cell = _cell_from_header(fr[1])
                    if cell is None:
                        print(f"[WARN] no Lattice in a frame of {f}; skip")
                        continue
                    sp, pos = _frames_to_species_pos(fr)
                    conv.append(("tuple", (cell, sp, pos)))
            except Exception as e:
                print(f"[WARN] xyz parse failed {f}: {e}")
                continue
        for k, item in enumerate(conv, start=1):
            dst = os.path.join(out, f"POSCAR_{next_idx}.vasp")
            if item[0] == "ase":
                _ase_write(dst, item[1], format="vasp", direct=False)
            else:
                cell, sp, pos = item[1]
                _write_poscar_vasp(dst, cell, sp, pos)
            print(f"[ABACUS] {base} frame {k} -> {os.path.basename(dst)}")
            next_idx += 1
    for f in sorted(glob.glob(os.path.join(src, "*.data"))):
        if not HAVE_ASE:
            print(f"[WARN] lammps-data {f} needs ase (use nepkit/HQmdstemkit env)")
            continue
        try:
            atoms = _ase_read(f, format="lammps-data")
        except Exception as e:
            print(f"[WARN] lammps-data read failed {f}: {e}")
            continue
        dst = os.path.join(out, f"POSCAR_{next_idx}.vasp")
        _ase_write(dst, atoms, format="vasp", direct=False)
        print(f"[ABACUS] {os.path.basename(f)} -> {os.path.basename(dst)}")
        next_idx += 1
    converted = sorted(glob.glob(os.path.join(out, "POSCAR_*.vasp")))
    if converted:
        names = ", ".join(os.path.basename(x) for x in converted)
        print(f"[ABACUS] POSCAR set ({len(converted)}): {names}")
    return out
def ask_path(prompt, default):
    if not sys.stdin.isatty():
        return default
    try:
        v = input(f"{prompt} [{default}]: ").strip()
    except EOFError:
        v = ""
    return v or default

def cmd_pretreat(args):
    if len(args) < 1:
        fail(USAGE)
    src = args[0]
    if os.path.isfile(src):
        src = os.path.dirname(os.path.abspath(src)) or "."
    prefix, out = "iter00", "."
    func = None
    data_name = None
    pp_dir = os.path.abspath(PP_DIR)
    orb_dir = os.path.abspath(ORB_DIR)
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--prefix" and i+1 < len(args): prefix = args[i+1]; i += 2
        elif a == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        elif a in ("--pp", "--pseudo-dir") and i+1 < len(args): pp_dir = args[i+1]; i += 2
        elif a in ("--orb", "--orbital-dir") and i+1 < len(args): orb_dir = args[i+1]; i += 2
        elif a == "--func" and i+1 < len(args): func = args[i+1].upper(); i += 2
        elif a == "--data-name" and i+1 < len(args): data_name = args[i+1]; i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    print(f"[ABACUS] bundled pseudo_dir : {os.path.abspath(PP_DIR)}")
    print(f"[ABACUS] bundled orbital_dir: {os.path.abspath(ORB_DIR)}")
    if func is None:
        func = ask_path("ABACUS functional (LDA/PBE/PBE-D3)", "PBE-D3").upper()
    if func not in ("LDA", "PBE", "PBE-D3"):
        fail(f"--func must be LDA, PBE or PBE-D3\n\n{USAGE}")
    dft_func = "PBE" if func in ("PBE", "PBE-D3") else "LDA"
    vdw = "vdw_method   d3_bj\n" if func == "PBE-D3" else ""
    print(f"[ABACUS] functional: {func}")
    os.makedirs(out, exist_ok=True)
    src = _convert_to_poscars(src, out)
    poscars = sorted(glob.glob(os.path.join(src, "POSCAR_*.vasp")))
    if not poscars:
        poscars = sorted(glob.glob(os.path.join(src, "*.vasp")))
    if not poscars:
        fail(f"no POSCAR_*.vasp in {src}\n\n{USAGE}")
    frames = []
    species_all = set()
    for f in poscars:
        lat, species, pos = read_poscar(f)
        frames.append((os.path.basename(f), lat, species, pos))
        species_all.update(species)
    os.makedirs(out, exist_ok=True)
    use_custom = pp_dir != os.path.abspath(PP_DIR) or orb_dir != os.path.abspath(ORB_DIR)
    if use_custom:
        pp_out = os.path.abspath(pp_dir)
        orb_out = os.path.abspath(orb_dir)
        print("[ABACUS] using custom pseudo/orbital dirs (no file copy)")
    else:
        data_name = data_name or ask_path("DFT data folder name under $HOME", "CuZn")
        data_dir = os.path.join(os.path.expanduser("~"), data_name.strip() or "CuZn")
        os.makedirs(data_dir, exist_ok=True)
        n = 0
        for el in sorted(species_all):
            pp_cands = sorted(glob.glob(os.path.join(PP_DIR, el + "_*.upf")) +
                              glob.glob(os.path.join(PP_DIR, el + "*.UPF")) +
                              glob.glob(os.path.join(PP_DIR, el + ".upf")))
            orb_cands = sorted(glob.glob(os.path.join(ORB_DIR, el + "_*.orb")))
            if pp_cands:
                shutil.copyfile(pp_cands[0], os.path.join(data_dir, os.path.basename(pp_cands[0])))
                n += 1
            if orb_cands:
                shutil.copyfile(orb_cands[0], os.path.join(data_dir, os.path.basename(orb_cands[0])))
                n += 1
        pp_out = data_dir
        orb_out = data_dir
        print(f"[ABACUS] data folder: {data_dir}")
        print(f"[ABACUS] copied {n} element files ({','.join(sorted(species_all))}) into {data_dir}")
    has_abacustest = shutil.which("abacustest") is not None
    for idx0, (name, lat, species, pos) in enumerate(frames, start=1):
        m = re.search(r"POSCAR_?(\d+)", name)
        idx = m.group(1) if m else str(idx0)
        d = os.path.join(out, f"{prefix}_{idx}")
        stru_done = False
        if has_abacustest:
            subprocess.run(["abacustest", "model", "inputs", "-f", os.path.join(src, name),
                            "--ftype", "poscar", "--lcao", "--jtype", "scf"],
                           cwd=src, check=False)
            ab0 = os.path.join(src, "000000")
            if os.path.isdir(ab0):
                if os.path.isdir(d):
                    shutil.rmtree(d)
                shutil.move(ab0, d)
                stru_path = os.path.join(d, "STRU")
                if os.path.isfile(stru_path):
                    lines = open(stru_path, encoding="utf-8-sig", errors="replace").read().splitlines()
                    idx_lat = next((i for i, l in enumerate(lines) if l.strip().startswith("LATTICE_CONSTANT")), len(lines))
                    with open(stru_path, "w", encoding="utf-8") as fh:
                        fh.write(_stru_header(species))
                        fh.write("\n".join(lines[idx_lat:]) + "\n")
                    stru_done = True
        if not stru_done:
            if not os.path.isdir(d):
                os.makedirs(d, exist_ok=True)
            write_stru(os.path.join(d, "STRU"), lat, species, pos)
        si = os.path.join(d, "struinfo.txt")
        if not os.path.isfile(si):
            with open(si, "w", encoding="utf-8") as fh:
                fh.write(name + "\n")
        with open(os.path.join(d, "INPUT"), "w", encoding="utf-8") as fh:
            fh.write(ABACUS_INPUT.format(pseudo=pp_out, orbital=orb_out, func=dft_func, vdw=vdw))
        print(f"{name} -> {d} (STRU/INPUT)  func={func}")
    sbatch_path = os.path.join(out, "run_array.sbatch2")
    with open(sbatch_path, "w", encoding="utf-8") as fh:
        fh.write(ABACUS_SBATCH.format(n=len(frames)))
    os.chmod(sbatch_path, 0o755)
    print(f"submission script: {sbatch_path}  (submit from {out})")

def cmd_extract(args):
    base = args[0] if args else "."
    out_xyz = "total_train.xyz"
    try:
        import dpdata
        from ase.io import write as ase_write
    except Exception as e:
        print(f"[WARN] dpdata/ase not available ({e}); fallback to batch_convert.sh")
        sh = os.path.join(EX_DIR, "batch_convert.sh")
        if os.path.isfile(sh):
            subprocess.run(["bash", sh], cwd=base, check=False)
        return
    dirs = []
    for root, subdirs, _ in os.walk(base):
        if "OUT.ABACUS" in subdirs:
            dirs.append(os.path.join(root, "OUT.ABACUS"))
    dirs.sort()
    if not dirs:
        print(f"no OUT.ABACUS found under {base}")
        return
    n = 0
    with open(out_xyz, "w", encoding="utf-8") as fh:
        pass
    for d in dirs:
        try:
            sys = dpdata.LabeledSystem(os.path.dirname(d), fmt="abacus/scf")
            traj = sys.to_ase_structure()
            atoms_list = traj if isinstance(traj, list) else [traj]
            for atoms in atoms_list:
                ase_write(out_xyz, atoms, format="extxyz", append=True)
                n += 1
        except Exception as e:
            print(f"[WARN] skip {d}: {e}")
    print(f"saved -> {out_xyz}  ({n} frames, NEP-ready extxyz)")

def cmd_shift(args):
    if len(args) < 2:
        fail(USAGE)
    conv = os.path.join(EX_DIR, "convert1ntest.py")
    if os.path.isfile(conv):
        subprocess.run([sys.executable, conv], cwd=os.path.dirname(os.path.abspath(args[0])), check=False)
    else:
        fail(f"convert1ntest.py not found\n\n{USAGE}")

def cmd_sbatch(args):
    print(ABACUS_SBATCH.format(n="N"))

def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "pretreat": cmd_pretreat(args)
    elif cmd == "extract": cmd_extract(args)
    elif cmd == "shift": cmd_shift(args)
    elif cmd == "sbatch": cmd_sbatch(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()