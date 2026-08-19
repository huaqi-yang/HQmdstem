#!/usr/bin/env python3
"""HQmdstemkit: CP2K single-point batch pretreatment / energy extraction.

Usage:
  hq_cp2k.py pretreat SRC [--func LDA|PBE|PBE-D3] [--elements Cu,Zn] [--prefix iter00]
                        [--out OUT] [--data DIR] [--basis FILE] [--potential FILE]
                        [--dftd3 FILE] [--queue gpu|cpu]
  hq_cp2k.py extract DIR [--csv cp2k_energies.csv]
  hq_cp2k.py sbatch [--queue gpu|cpu]
"""
import os, sys, glob, re, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_batch import read_poscar
from hq_common import fail, hq_home, read_extxyz_frames

USAGE = __doc__
CP2K_EX = os.path.join(hq_home(), "examples", "cp2k")
DATA_DIR = os.path.join(CP2K_EX, "data")

# CP2K GTH potential/basis q-values for common elements.
ELEM_Q = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 4, "N": 5, "O": 6,
          "F": 7, "Ne": 8, "Na": 9, "Mg": 10, "Al": 11, "Si": 4, "P": 5,
          "S": 6, "Cl": 7, "Ar": 8, "K": 9, "Ca": 10, "Ti": 12, "V": 13,
          "Cr": 14, "Mn": 15, "Fe": 16, "Co": 17, "Ni": 18, "Cu": 11, "Zn": 12,
          "Ge": 4, "Se": 6, "Sr": 10, "Zr": 12, "Mo": 14, "Ag": 11, "Sn": 4,
          "Te": 6, "Ba": 10, "Hf": 12, "W": 14, "Pt": 10, "Au": 11, "Pb": 4}
FUNC_INFO = {
    "LDA":   {"xc": "      &XC_FUNCTIONAL\n        &LDA\n          FUNCTIONAL PADE\n        &END LDA\n      &END XC_FUNCTIONAL\n",
              "pot": "GTH-LDA", "label": "LDA"},
    "PBE":   {"xc": "      &XC_FUNCTIONAL\n        &PBE\n        &END PBE\n      &END XC_FUNCTIONAL\n",
              "pot": "GTH-PBE", "label": "PBE"},
    "PBE-D3": {"xc": "      &XC_FUNCTIONAL\n        &PBE\n        &END PBE\n      &END XC_FUNCTIONAL\n",
              "pot": "GTH-PBE", "label": "PBE-D3"},
}

GPU_ARRAY = """#!/bin/bash
#SBATCH --job-name=CuZn-cp2k
#SBATCH --partition=4V100
#SBATCH --nodes=1
#SBATCH --ntasks=4          # 1 node x 4 GPUs
#SBATCH --gpus-per-node=4
#SBATCH --qos=flood-1o2gpu
#SBATCH --array=1-{n}%50

# Do not modify CUDA-MPS / rank-map settings unless you know what you are doing.
source /opt/sai_config/mps_mapping.d/${{SLURM_JOB_PARTITION}}.bash
export OMP_NUM_THREADS=$((CORES_PER_GPU/RANKS_PER_GPU))
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load cp2k/2025.1-cuda12.4-sm70-auto
WORK_DIR="$PWD/iter00_${{SLURM_ARRAY_TASK_ID}}"
cd $WORK_DIR
mpirun -np $SLURM_NTASKS --map-by $MAP_OPT cp2k.psmp -i input.inp -o output.log
"""

CPU_ARRAY = """#!/bin/bash
#SBATCH --job-name=CuZn-cp2k
#SBATCH --partition=D9654
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --qos=huge-cpu
#SBATCH --array=1-{n}%50

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load cp2k/2025.1-cpu-auto
WORK_DIR="$PWD/iter00_${{SLURM_ARRAY_TASK_ID}}"
cd $WORK_DIR
mpirun -np $SLURM_NTASKS --map-by ppr:$((8/OMP_NUM_THREADS)):l3cache:pe=$OMP_NUM_THREADS -mca coll_hcoll_enable 0 cp2k.psmp -i input.inp -o output.log
"""

DATA_URLS = {
    "BASIS_MOLOPT_UZH": "https://github.com/cp2k/cp2k/raw/master/data/BASIS_MOLOPT_UZH",
    "POTENTIAL_UZH":    "https://github.com/cp2k/cp2k/raw/master/data/POTENTIAL_UZH",
    "dftd3.dat":        "https://github.com/cp2k/cp2k/raw/master/data/dftd3.dat",
}


def _cell_from_frames(frames):
    """Extract 3x3 lattice from an extxyz header, return None if absent."""
    hdr = frames[1]
    m = re.search(r"Lattice\s*=\s*[\"']?([-+0-9.eE\s]+)[\"']?", hdr)
    if m:
        vals = m.group(1).split()
        if len(vals) >= 9:
            return np.array([[float(x) for x in vals[i*3:i*3+3]] for i in range(3)])
    m = re.search(r"lattice\s*=\s*[\"']?([-+0-9.eE\s]+)[\"']?", hdr)
    if m:
        vals = m.group(1).split()
        if len(vals) >= 9:
            return np.array([[float(x) for x in vals[i*3:i*3+3]] for i in range(3)])
    return None


def _frames_to_xyz(frames):
    """Convert extxyz frame to (species_list, pos_array) without lattice."""
    n = frames[0]
    species, coords = [], []
    for al in frames[2:2 + n]:
        c = al.split()
        if len(c) >= 4:
            species.append(c[0].capitalize())
            coords.append([float(c[1]), float(c[2]), float(c[3])])
    return species, np.array(coords)


def _data_file(explicit, data_dir, name):
    if explicit and os.path.isfile(explicit):
        return os.path.abspath(explicit)
    p = os.path.join(data_dir, name)
    if os.path.isfile(p):
        return p
    return None


def ask_choice(prompt, default):
    if not sys.stdin.isatty():
        return default
    try:
        v = input(f"{prompt} [{default}]: ").strip()
    except EOFError:
        v = ""
    return v or default

def _build_input(func, cell, species, basis_path, pot_path, dftd3_path, basis_name, pot_prefix):
    kinds = []
    for el in sorted(set(species)):
        q = ELEM_Q.get(el)
        basis = f"{basis_name}-q{q}" if q else basis_name
        pot = f"{pot_prefix}-q{q}" if q else pot_prefix
        kinds.append(f"    &KIND {el}\n      ELEMENT {el}\n"
                     f"      BASIS_SET {basis}\n      POTENTIAL {pot}\n    &END KIND\n")
    basis_line = basis_path if basis_path else "./BASIS_MOLOPT_UZH"
    pot_line = pot_path if pot_path else "./POTENTIAL_UZH"
    dftd3_line = dftd3_path if dftd3_path else "./dftd3.dat"
    vdw = ""
    if func == "PBE-D3":
        vdw = ("      &VDW_POTENTIAL\n        POTENTIAL_TYPE PAIR_POTENTIAL\n"
               "        &PAIR_POTENTIAL\n"
               f"          PARAMETER_FILE_NAME {dftd3_line}\n"
               "          TYPE DFTD3\n          DAMPING_TYPE BJ\n          REFERENCE_FUNCTIONAL PBE\n"
               "        &END PAIR_POTENTIAL\n      &END VDW_POTENTIAL\n")
    a, b, c = cell
    txt = f"""&GLOBAL
  PROJECT cp2k
  PRINT_LEVEL MEDIUM
  RUN_TYPE ENERGY_FORCE
  EXTENDED_FFT_LENGTHS T
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &CELL
      A     {a[0]:14.8f} {a[1]:14.8f} {a[2]:14.8f}
      B     {b[0]:14.8f} {b[1]:14.8f} {b[2]:14.8f}
      C     {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME pos.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
{''.join(kinds)}  &END SUBSYS

  &DFT
    BASIS_SET_FILE_NAME  {basis_line}
    POTENTIAL_FILE_NAME  {pot_line}
    CHARGE 0
    MULTIPLICITY 1
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &POISSON
      PERIODIC XYZ
      PSOLVER PERIODIC
    &END POISSON
    &XC
{FUNC_INFO[func]['xc']}{vdw}    &END XC
    &MGRID
      CUTOFF 500
      REL_CUTOFF 60
    &END MGRID
    &SCF
      SCF_GUESS ATOMIC
      MAX_SCF 100
      EPS_SCF 1.0E-07
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
    &END SCF
  &END DFT
  &PRINT
    &FORCES ON
    &END FORCES
    &STRESS ON
    &END STRESS
  &END PRINT
&END FORCE_EVAL
"""
    return txt


def _write_pos_xyz(path, species, pos):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{len(species)}\ncp2k single point\n")
        for s, p in zip(species, pos):
            f.write(f"{s:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")


def cmd_pretreat(args):
    if len(args) < 1:
        fail(USAGE)
    src = args[0]
    if os.path.isfile(src):
        src = os.path.dirname(os.path.abspath(src)) or "."
    func = None
    elements = None
    data_name = None
    prefix, out = "iter00", "."
    data_dir = DATA_DIR
    basis_file = pot_file = dftd3_file = None
    queue = "gpu"
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--func" and i+1 < len(args): func = args[i+1].upper(); i += 2
        elif a == "--prefix" and i+1 < len(args): prefix = args[i+1]; i += 2
        elif a == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        elif a == "--data" and i+1 < len(args): data_dir = args[i+1]; i += 2
        elif a == "--basis" and i+1 < len(args): basis_file = args[i+1]; i += 2
        elif a == "--potential" and i+1 < len(args): pot_file = args[i+1]; i += 2
        elif a == "--dftd3" and i+1 < len(args): dftd3_file = args[i+1]; i += 2
        elif a == "--queue" and i+1 < len(args): queue = args[i+1].lower(); i += 2
        elif a == "--elements" and i+1 < len(args): elements = args[i+1]; i += 2
        elif a == "--data-name" and i+1 < len(args): data_name = args[i+1]; i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    basis_path = _data_file(basis_file, data_dir, "BASIS_MOLOPT_UZH")
    pot_path = _data_file(pot_file, data_dir, "POTENTIAL_UZH")
    d3_path = _data_file(dftd3_file, data_dir, "dftd3.dat")
    missing = [n for n, p in [("BASIS_MOLOPT_UZH", basis_path), ("POTENTIAL_UZH", pot_path), ("dftd3.dat", d3_path)] if not p]
    os.makedirs(out, exist_ok=True)
    if missing:
        print(f"[CP2K] data dir used: {data_dir}")
        print("[CP2K] missing data files: " + ", ".join(missing))
        print("[CP2K] download them into examples/cp2k/data/ :")
        for n in missing:
            print(f"  {DATA_URLS[n]}")
    else:
        data_name = data_name or ask_choice("DFT data folder name under $HOME", "CuZn")
        home_data = os.path.join(os.path.expanduser("~"), data_name.strip() or "CuZn")
        os.makedirs(home_data, exist_ok=True)
        for spath, name in [(basis_path, "BASIS_MOLOPT_UZH"), (pot_path, "POTENTIAL_UZH"), (d3_path, "dftd3.dat")]:
            if spath and os.path.isfile(spath):
                shutil.copyfile(spath, os.path.join(home_data, name))
        basis_path = os.path.join(home_data, "BASIS_MOLOPT_UZH")
        pot_path = os.path.join(home_data, "POTENTIAL_UZH")
        d3_path = os.path.join(home_data, "dftd3.dat")
        print(f"[CP2K] data folder: {home_data}")
        print(f"[CP2K] copied BASIS_MOLOPT_UZH / POTENTIAL_UZH / dftd3.dat into {home_data}")
    # collect structures
    frames = []
    poscars = sorted(glob.glob(os.path.join(src, "POSCAR_*.vasp")))
    if not poscars:
        poscars = sorted(glob.glob(os.path.join(src, "*.vasp")))
    xyzs = sorted(glob.glob(os.path.join(src, "*.xyz")))
    if poscars:
        for f in poscars:
            lat, species, pos = read_poscar(f)
            frames.append((os.path.basename(f), lat, species, pos))
    elif xyzs:
        for xf in xyzs:
            all_frames = read_extxyz_frames(xf)
            if not all_frames:
                continue
            for fr in all_frames:
                cell = _cell_from_frames(fr)
                if cell is None:
                    fail(f"no Lattice=... in {xf}; CP2K needs a periodic cell")
                species, pos = _frames_to_xyz(fr)
                frames.append((os.path.basename(xf), cell, species, pos))
    if not frames:
        fail(f"no POSCAR_*.vasp or *.xyz frames found in {src}\n\n{USAGE}")
    all_species = sorted({s for _, _, species, _ in frames for s in species})
    if func is None:
        func = ask_choice("Functional (LDA/PBE/PBE-D3)", "PBE-D3").upper()
    if func not in FUNC_INFO:
        fail(f"--func must be LDA, PBE or PBE-D3\n\n{USAGE}")
    if elements is None:
        elements = ask_choice("Element species (comma separated, auto)", ",".join(all_species))
    sel = set(e.capitalize() for e in re.split(r"[,，\s]+", elements) if e) if elements else set(all_species)
    if not sel:
        sel = set(all_species)
    print(f"[CP2K] functional: {FUNC_INFO[func]['label']}   elements: {','.join(sorted(sel))}")
    os.makedirs(out, exist_ok=True)
    basis_name = "DZVP-MOLOPT-SR-GTH"
    pot_prefix = FUNC_INFO[func]["pot"]
    for idx, (name, cell, species, pos) in enumerate(frames, start=1):
        keep = [s in sel for s in species]
        sp2 = [s for s, k in zip(species, keep) if k]
        pos2 = np.array([p for p, k in zip(pos, keep) if k])
        if len(sp2) != len(species):
            print(f"[WARN] {name}: dropping {len(species)-len(sp2)} atoms outside element selection {sorted(sel)}")
        d = os.path.join(out, f"{prefix}_{idx}")
        os.makedirs(d, exist_ok=True)
        _write_pos_xyz(os.path.join(d, "pos.xyz"), sp2, pos2)
        inp = _build_input(func, cell, sp2, basis_path, pot_path, d3_path,
                           basis_name, pot_prefix)
        with open(os.path.join(d, "input.inp"), "w", encoding="utf-8") as f:
            f.write(inp)
        print(f"{name} -> {d} (pos.xyz/input.inp)  func={FUNC_INFO[func]['label']}")
    n = len(frames)
    gpu_sbatch = os.path.join(out, "run_cp2k_gpu_array.sbatch")
    cpu_sbatch = os.path.join(out, "run_cp2k_cpu_array.sbatch")
    with open(gpu_sbatch, "w", encoding="utf-8") as f:
        f.write(GPU_ARRAY.format(n=n))
    with open(cpu_sbatch, "w", encoding="utf-8") as f:
        f.write(CPU_ARRAY.format(n=n))
    os.chmod(gpu_sbatch, 0o755)
    os.chmod(cpu_sbatch, 0o755)
    print(f"submission scripts: {gpu_sbatch} (GPU), {cpu_sbatch} (CPU)")
    print(f"submit on SAI: sbatch run_cp2k_gpu_array.sbatch   (or run_cp2k_cpu_array.sbatch)")
    if missing:
        print("[CP2K] NOTE: download the missing data files first (see URLs above).")


def _read_stress(log):
    try:
        lines = open(log, encoding="utf-8-sig", errors="replace").read().splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        if "Analytical stress tensor" in line:
            for ln in lines[i+1:i+4]:
                nums = []
                for tok in ln.replace(",", " ").split():
                    try:
                        nums.append(float(tok))
                    except ValueError:
                        pass
                if len(nums) >= 6:
                    return nums[:6]
    return None


def _cell_vectors(d):
    inp = os.path.join(d, "input.inp")
    if not os.path.isfile(inp):
        return None
    vecs = []
    for line in open(inp, encoding="utf-8-sig", errors="replace"):
        s = line.split()
        if len(s) == 4 and s[0] in ("A", "B", "C"):
            try:
                vecs.append([float(x) for x in s[1:]])
            except ValueError:
                pass
    return np.array(vecs) if len(vecs) == 3 else None


def _cell_volume(d):
    vecs = _cell_vectors(d)
    if vecs is None:
        return None
    return abs(float(np.linalg.det(vecs)))


def cmd_extract(args):
    base = args[0] if args else "."
    csv_out = "cp2k_energies.csv"
    if "--csv" in args:
        csv_out = args[args.index("--csv") + 1]
    rows = []
    for d in sorted(glob.glob(os.path.join(base, "iter00_*"))):
        log = os.path.join(d, "output.log")
        if not os.path.isfile(log):
            continue
        energy = None
        for line in open(log, encoding="utf-8-sig", errors="replace"):
            if "Total FORCE_EVAL ( QS ) energy" in line:
                parts = line.split()
                for tok in reversed(parts):
                    try:
                        energy = float(tok)
                        break
                    except ValueError:
                        continue
            if energy is not None:
                break
        stress = _read_stress(log)
        V = _cell_volume(d)
        rows.append((os.path.basename(d), energy, stress, V, _cell_vectors(d)))
    if not rows:
        print(f"no energies found under {base} (look for iter00_*/output.log)")
        return
    with open(csv_out, "w", encoding="utf-8") as f:
        f.write("iter,energy_au,energy_eV,"
                "stress_xx_GPa,stress_yy_GPa,stress_zz_GPa,stress_xy_GPa,stress_xz_GPa,stress_yz_GPa,"
                "virial_xx_eV,virial_yy_eV,virial_zz_eV,virial_xy_eV,virial_xz_eV,virial_yz_eV\n")
        for name, energy, stress, V, _cell in rows:
            e_au = f"{energy:.12f}" if energy is not None else ""
            e_ev = f"{energy * 27.211386245988:.8f}" if energy is not None else ""
            if stress:
                s6 = [f"{x:.8f}" for x in stress]
                if V:
                    v6 = [f"{-x * V * 0.006241509:.8f}" for x in stress]
                else:
                    v6 = [""] * 6
            else:
                s6 = [""] * 6
                v6 = [""] * 6
            f.write(f"{name},{e_au},{e_ev}," + ",".join(s6 + v6) + "\n")
    print(f"saved -> {csv_out}  ({len(rows)} energies, stress + virial columns)")
    xyz_out = "cp2k_train.xyz"
    n_frames = 0
    with open(xyz_out, "w", encoding="utf-8") as f:
        for name, energy, stress, V, cell in rows:
            pos_file = os.path.join(base, name, "pos.xyz")
            if energy is None or not os.path.isfile(pos_file):
                continue
            lines = open(pos_file, encoding="utf-8-sig", errors="replace").read().splitlines()
            try:
                n = int(lines[0].split()[0])
            except Exception:
                continue
            atoms = lines[2:2 + n]
            e_ev = energy * 27.211386245988
            if cell is not None:
                lat = " ".join(f"{x:.10f}" for row in cell for x in row)
                hdr = f'Lattice="{lat}" Properties=species:S:1:pos:R:3 energy={e_ev:.8f} pbc="T T T"'
            else:
                hdr = f'Properties=species:S:1:pos:R:3 energy={e_ev:.8f} pbc="T T T"'
            if stress and V:
                vir = [-x * V * 0.006241509 for x in stress]
                hdr += ' virial="' + " ".join(f"{x:.8f}" for x in vir) + '"'
            f.write(str(n) + "\n" + hdr + "\n")
            for al in atoms:
                f.write(al + "\n")
            n_frames += 1
    print(f"saved -> {xyz_out}  ({n_frames} NEP-ready extxyz frames)")


def cmd_sbatch(args):
    queue = "gpu"
    if "--queue" in args:
        queue = args[args.index("--queue") + 1]
    if queue == "cpu":
        print(CPU_ARRAY.format(n="N"))
    else:
        print(GPU_ARRAY.format(n="N"))


def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "pretreat": cmd_pretreat(args)
    elif cmd == "extract": cmd_extract(args)
    elif cmd == "sbatch": cmd_sbatch(args)
    else: fail(USAGE)


if __name__ == "__main__":
    main()