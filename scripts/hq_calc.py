#!/usr/bin/env python3
"""HQmdstemkit: calculators (cohesive/EOS, shear, stacking-fault).

Usage:
  hq_calc.py cohesive DIR [--elements Cu,Zn] [--min 0.9 --max 1.2 --n 6] [--out eos.png]
  hq_calc.py shear DIR            # TODO tomorrow
  hq_calc.py stacking-fault DIR   # TODO tomorrow
"""
import os, sys, re, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, hq_home

USAGE = __doc__
CALC_DIR = os.path.join(hq_home(), "examples", "calculators")


def ask_choice(prompt, default):
    if not sys.stdin.isatty():
        return default
    try:
        v = input(f"{prompt} [{default}]: ").strip()
    except EOFError:
        v = ""
    return v or default


def _copy_potential(d, elements):
    sel = set(e.capitalize() for e in re.split(r"[,，\s]+", elements) if e)
    if {"Cu", "Zn"} <= sel:
        src = os.path.join(hq_home(), "examples", "gpumd", "nep.txt")
        label = "CuZn dedicated"
    else:
        src = os.path.join(hq_home(), "examples", "nep", "nep89_20250409.txt")
        label = "generic NEP89"
    if os.path.isfile(src):
        shutil.copyfile(src, os.path.join(d, "nep.txt"))
    print(f"[CALC] potential ({label}) -> {os.path.join(d, 'nep.txt')}")


def _model_volume(d):
    mxyz = os.path.join(d, "model.xyz")
    if not os.path.isfile(mxyz):
        return None
    with open(mxyz, encoding="utf-8-sig", errors="replace") as f:
        lines = f.readlines()
    if not lines:
        return None
    hdr = lines[1] if len(lines) > 1 and lines[0].strip().lstrip("-").isdigit() else lines[0]
    m = re.search(r"Lattice\s*=\s*[\"']?([-+0-9.eE\s]+)[\"']?", hdr)
    if not m:
        return None
    vals = [float(x) for x in m.group(1).split()]
    if len(vals) < 9:
        return None
    lat = np.array([vals[i*3:i*3+3] for i in range(3)])
    return abs(float(np.linalg.det(lat)))


def cmd_cohesive(args):
    d = args[0] if args else "."
    elements = None; smin, smax, n = 0.9, 1.2, 6
    out = "energy_vs_volume_cohesive.png"
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--elements" and i+1 < len(args): elements = args[i+1]; i += 2
        elif a == "--min" and i+1 < len(args): smin = float(args[i+1]); i += 2
        elif a == "--max" and i+1 < len(args): smax = float(args[i+1]); i += 2
        elif a == "--n" and i+1 < len(args): n = int(args[i+1]); i += 2
        elif a == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    os.makedirs(d, exist_ok=True)
    if elements is None:
        elements = ask_choice("Element species (comma separated)", "Cu,Zn")
    run_in = os.path.join(d, "run.in")
    if not os.path.isfile(run_in):
        with open(run_in, "w") as f:
            f.write(f"potential nep.txt\ncompute_cohesive {smin} {smax} {n}\n")
        print(f"[CALC] {run_in} written")
    if not os.path.isfile(os.path.join(d, "model.xyz")):
        shutil.copyfile(os.path.join(hq_home(), "examples", "gpumd", "model.xyz"),
                        os.path.join(d, "model.xyz"))
    _copy_potential(d, elements)
    coh = os.path.join(d, "cohesive.out")
    if not os.path.isfile(coh):
        print("[CALC] run.in/nep.txt/model.xyz ready; next: gpumd in this dir, then rerun 4 408")
        return
    data = np.loadtxt(coh)
    scale, energy = data[:, 0], data[:, 1]
    V0 = _model_volume(d)
    if V0 is None:
        print("[CALC] model.xyz missing/without Lattice; cannot convert scale -> volume")
        return
    vol = V0 * scale**3
    order = np.argsort(np.abs(scale - 1.0))
    k = max(5, len(scale) // 2)
    idx = np.sort(order[:k])
    if len(idx) >= 3:
        coeff = np.polyfit(vol[idx], energy[idx], 2)
        vmin = -coeff[1] / (2 * coeff[0])
        emin = np.polyval(coeff, vmin)
        smin_fit = (vmin / V0) ** (1.0 / 3.0)
        print(f"[CALC] E_min = {emin:.6f} eV")
        print(f"[CALC] V_min = {vmin:.4f} A^3")
        print(f"[CALC] scale_min = {smin_fit:.5f}  (cubic a = {vmin**(1.0/3.0):.5f} A)")
    else:
        print("[CALC] too few points for parabola fit")
    np.savetxt(os.path.join(d, "energy_vs_volume.txt"),
               np.column_stack((vol, energy)), header="Volume(A^3) Energy(eV)")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
        ax.plot(vol, energy, "o-", lw=2.0, ms=4, color="#1f77b4")
        if len(idx) >= 3:
            vv = np.linspace(vol.min(), vol.max(), 200)
            ax.plot(vv, np.polyval(coeff, vv), "--", color="#ff7f0e", lw=1.5)
            ax.scatter([vmin], [emin], s=60, marker="*", color="red", zorder=5)
        ax.set_xlabel("Volume ($\\AA^3$)")
        ax.set_ylabel("Energy (eV)")
        ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout()
        fig.savefig(os.path.join(d, out), dpi=300)
        print(f"[CALC] plot -> {os.path.join(d, out)}")
    except Exception as e:
        print(f"[WARN] plotting skipped: {e}")


def cmd_todo(name, args):
    d = args[0] if args else "."
    print(f"[CALC] {name}: not implemented yet (tomorrow). Directory: {d}")


def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "cohesive": cmd_cohesive(args)
    elif cmd == "shear": cmd_todo("Shear", args)
    elif cmd == "stacking-fault": cmd_todo("Stacking fault", args)
    else: fail(USAGE)


if __name__ == "__main__":
    main()