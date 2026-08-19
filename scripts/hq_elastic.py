#!/usr/bin/env python3
"""HQmdstemkit: elastic constants (strain fluctuation + 0K compute_elastic).

Usage:
  hq_elastic.py strain DIR [--T 300] [--skip 1000] [--slices 10]
  hq_elastic.py 0k DIR [--strain 0.01] [--run]
  hq_elastic.py plot CSV [--out elastic.png]
  hq_elastic.py born ELASTIC_OUT
"""
import glob, os, shutil, sys, math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, parse_thermo_out, load_config, hq_home

USAGE = __doc__

def _angle(v1, v2):
    return math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def analyze_strain_fluctuation(thermo, T, slices):
    """Cubic elastic constants C11,C12,C44 from box-vector fluctuations."""
    kB = 1.38064852e-23
    n = len(thermo["ax"])
    ax = np.asarray(thermo["ax"], dtype=float); ay = np.asarray(thermo["ay"], dtype=float)
    az = np.asarray(thermo["az"], dtype=float); bx = np.asarray(thermo["bx"], dtype=float)
    by = np.asarray(thermo["by"], dtype=float); bz = np.asarray(thermo["bz"], dtype=float)
    cx = np.asarray(thermo["cx"], dtype=float); cy = np.asarray(thermo["cy"], dtype=float)
    cz = np.asarray(thermo["cz"], dtype=float)
    V = np.mean(np.array(ax) * np.array(by) * np.array(cz))  # A^3
    scale = 100.0 / (T * kB) * V  # GPa^-1
    per = n // slices
    C = []
    for i in range(slices):
        s = slice(i * per, (i + 1) * per)
        aa = ax[s]/np.mean(ax[s]) - 1
        bb = by[s]/np.mean(by[s]) - 1
        cc = cz[s]/np.mean(cz[s]) - 1
        alpha = np.array([_angle([bx[j], by[j], bz[j]], [cx[j], cy[j], cz[j]]) for j in range(s.start, s.stop)])
        beta  = np.array([_angle([ax[j], ay[j], az[j]], [cx[j], cy[j], cz[j]]) for j in range(s.start, s.stop)])
        gamma = np.array([_angle([ax[j], ay[j], az[j]], [bx[j], by[j], bz[j]]) for j in range(s.start, s.stop)])
        e11, e22, e33 = aa, bb, cc
        e23 = (alpha - math.pi/2)/2
        e13 = (beta  - math.pi/2)/2
        e12 = (gamma - math.pi/2)/2
        cov = lambda x, y: np.mean(x*y) - np.mean(x)*np.mean(y)
        S11 = (cov(e11,e11)+cov(e22,e22)+cov(e33,e33))/3
        S12 = (cov(e11,e22)+cov(e11,e33)+cov(e22,e33))/3
        S44 = 4*(cov(e23,e23)+cov(e13,e13)+cov(e12,e12))/3
        S = np.array([[S11,S12,S12,0,0,0],[S12,S11,S12,0,0,0],[S12,S12,S11,0,0,0],
                      [0,0,0,S44,0,0],[0,0,0,0,S44,0],[0,0,0,0,0,S44]]) * scale
        try:
            Cpq = np.linalg.inv(S)
        except np.linalg.LinAlgError:
            continue
        C.append([Cpq[1,1], Cpq[1,2], Cpq[4,4]])
    C = np.array(C)
    if C.size == 0:
        raise RuntimeError("no valid slices: the compliance matrix was singular "
                           "(box/angles do not fluctuate enough)")
    return C.mean(axis=0), C.std(axis=0)

def cmd_strain(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    T = 300; skip = 1000; slices = 10
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--T" and i+1 < len(args): T = float(args[i+1]); i += 2
        elif a == "--skip" and i+1 < len(args): skip = int(args[i+1]); i += 2
        elif a == "--slices" and i+1 < len(args): slices = int(args[i+1]); i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    tp = os.path.join(d, "thermo.out")
    if not os.path.isfile(tp):
        fail(f"thermo.out not found in {d}\n\n{USAGE}")
    thermo = parse_thermo_out(tp, skip_frames=skip)
    try:
        avg, err = analyze_strain_fluctuation(thermo, T, slices)
    except RuntimeError as e:
        fail(f"{e}\n\n{USAGE}")
    print(f"T={T:.0f} K  strain-fluctuation elastic constants (GPa):")
    print(f"  C11 = {avg[0]:8.2f} +/- {err[0]:.2f}")
    print(f"  C12 = {avg[1]:8.2f} +/- {err[1]:.2f}")
    print(f"  C44 = {avg[2]:8.2f} +/- {err[2]:.2f}")

def cmd_0k(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    strain = 0.01; do_run = False
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--strain" and i+1 < len(args): strain = float(args[i+1]); i += 2
        elif a == "--run": do_run = True; i += 1
        else: fail(f"unknown option {a}\n\n{USAGE}")
    for f in ("model.xyz", "nep.txt"):
        if not os.path.isfile(os.path.join(d, f)):
            fail(f"missing {f} in {d}\n\n{USAGE}")
    with open(os.path.join(d, "run.in"), "w") as f:
        f.write(f"potential   ./nep.txt\n\nminimize fire 1.0e-5 1000 1\n"
                f"compute_elastic {strain:g}\n")
    print(f"0K compute_elastic run.in written in {d} (strain={strain:g})")
    if do_run:
        import subprocess
        subprocess.run(["gpumd"], cwd=d, check=True)
        print("GPUMD finished; see elastic.out")

def cmd_plot(args):
    if len(args) < 1:
        fail(USAGE)
    csv_path = args[0]
    out = "elastic_constants.png"
    if "--out" in args:
        out = args[args.index("--out") + 1]
    if not os.path.isfile(csv_path):
        fail(f"CSV not found: {csv_path}\n\nUSAGE: {USAGE}")
    import csv
    rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig", errors="replace")))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scales = sorted(set(float(r["Scale"]) for r in rows))
    cols = ["C11", "C22", "C33", "C44", "C55", "C66", "C12", "C13", "C23"]
    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
    for sc in scales:
        sub = [r for r in rows if float(r["Scale"]) == sc]
        T = [float(r["Temperature(K)"]) for r in sub]
        for col in cols:
            ax.plot(T, [float(r[col]) for r in sub], "o-", lw=1.2, ms=3,
                    label=f"Scale={sc:g} {col}")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Elastic constants (GPa)")
    ax.legend(fontsize=6, ncol=2, frameon=False)
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"saved -> {out}")

def cmd_born(args):
    if len(args) < 1:
        fail(USAGE)
    p = args[0]
    if not os.path.isfile(p):
        fail(f"elastic.out not found: {p}\n\n{USAGE}")
    C = np.loadtxt(p, comments="#")
    if C.shape != (6, 6):
        fail("elastic.out must be a 6x6 matrix\n\n" + USAGE)
    # Born stability for the general triclinic case: principal minors of C > 0
    ok = True
    for k in range(1, 7):
        if np.linalg.det(C[:k, :k]) <= 0:
            ok = False
            break
    print("C matrix (GPa):")
    print(C)
    print("Born stability (principal minors > 0):", "STABLE" if ok else "NOT STABLE")


BORN_COLORS = ["#4CAF50", "#8BC34A", "#FF9800", "#E53935"]
BORN_MARKERS = ["o", "s", "^", "D"]

def cmd_born_plot(args):
    """3x2 Born-stability plot from elastic_constants_raw_data.csv."""
    csv_path = None
    out = "Born_Stability_3x2.png"
    i = 0
    while i < len(args):
        if args[i] == "--csv" and i+1 < len(args): csv_path = args[i+1]; i += 2
        elif args[i] == "--out" and i+1 < len(args): out = args[i+1]; i += 2
        else: fail(f"unknown option {args[i]}\n\n{USAGE}")
    if csv_path is None:
        cfg_csv = load_config().get("elastic_raw_csv", "")
        cands = [cfg_csv] if cfg_csv else []
        cands += [os.path.join(hq_home(), "examples", "elastic_constants_raw_data.csv")]
        for cand in cands:
            if os.path.isfile(cand):
                csv_path = cand; break
    if not csv_path or not os.path.isfile(csv_path):
        fail(f"CSV not found; use --csv PATH\n\n{USAGE}")
    import csv
    rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig", errors="replace")))
    scales = sorted(set(float(r["Scale"]) for r in rows), key=float)
    temps = sorted(set(float(r["Temperature(K)"]) for r in rows))
    def series(scale, key):
        sub = [r for r in rows if float(r["Scale"]) == scale]
        y = []
        for t in temps:
            hit = [r for r in sub if float(r["Temperature(K)"]) == t]
            y.append(float(hit[0][key]) if hit else np.nan)
        return np.array(y)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(20/2.54, 22/2.54), sharex=True)
    fig.subplots_adjust(left=0.10, bottom=0.06, right=0.98, top=0.97, hspace=0.0, wspace=0.25)
    items = [("C11", r"$C_{11}$ (GPa)", 1e-2), ("C44", r"$C_{44}$ (GPa)", 1e-2),
             ("C55", r"$C_{55}$ (GPa)", 1e-3), ("C66", r"$C_{66}$ (GPa)", 1e-3)]
    for ax in axes.flat:
        ax.axhline(0, color="#999999", ls="--", lw=1.0)
    for idx, scale in enumerate(scales):
        color = BORN_COLORS[idx % len(BORN_COLORS)]
        marker = BORN_MARKERS[idx % len(BORN_MARKERS)]
        for j, (key, lab, sc) in enumerate(items):
            ax = axes.flat[j]
            y = series(scale, key) * sc
            m = ~np.isnan(y)
            ax.plot(np.array(temps)[m], y[m], color=color, marker=marker, ms=6,
                    markerfacecolor="none", lw=2.0, label=f"SCRAPS-{scale:g}" if j == 0 else None)
            ax.set_ylabel(lab)
        c11 = series(scale, "C11"); c22 = series(scale, "C22"); c33 = series(scale, "C33")
        c12 = series(scale, "C12"); c13 = series(scale, "C13"); c23 = series(scale, "C23")
        det2 = c11*c22 - c12**2
        det3 = c11*c22*c33 + 2*c12*c13*c23 - c11*c23**2 - c22*c13**2 - c33*c12**2
        for j, (y, sc, lab) in enumerate([(det2, 1e-5, r"$C_{11}C_{22}-C_{12}^2$ ($\times$10$^5$)"),
                                          (det3, 1e-7, r"Born $D_3$ ($\times$10$^7$)")]):
            ax = axes.flat[4 + j]
            m = ~np.isnan(y)
            ax.plot(np.array(temps)[m], y[m]*sc, color=color, marker=marker, ms=6,
                    markerfacecolor="none", lw=2.0, label=None)
            ax.set_ylabel(lab)
    for k, ax in enumerate(axes.flat):
        ax.tick_params(direction="in", top=True, right=True)
        ax.set_xlim(0, 700)
        if k < 4:
            ax.tick_params(axis="x", labelbottom=False)
        else:
            ax.set_xlabel("Temperature (K)")
    axes.flat[0].legend(fontsize=12, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    print(f"saved -> {out}")

def cmd_auto(args):
    """Auto elastic inputs: detect local model.xyz/nep.txt/run.in, else bundled example,
    then choose method 1) GPUMD 0K  2) GPUMD+NEP89  3) LAMMPS EAM."""
    d = None; method = None
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--method" and i+1 < len(args): method = args[i+1]; i += 2
        elif not a.startswith("-") and d is None: d = a; i += 1
        else: fail(f"unknown option {a}\n\n{USAGE}")
    d = d or os.getcwd()
    os.makedirs(d, exist_ok=True)
    src = load_config().get("elastic_example_dir") or os.path.join(hq_home(), "examples", "elastic")
    have = all(os.path.isfile(os.path.join(d, f)) for f in ("model.xyz", "nep.txt", "run.in"))
    if not have:
        for f in ("model.xyz", "nep.txt", "run.in"):
            s = os.path.join(src, f)
            if os.path.isfile(s) and not os.path.isfile(os.path.join(d, f)):
                shutil.copyfile(s, os.path.join(d, f))
        print(f"local elastic inputs missing -> copied bundled example from {src}")
    if method is None:
        try:
            method = input("Elastic method? 1) GPUMD 0K  2) GPUMD+NEP89  3) LAMMPS EAM  [1]: ").strip()
        except EOFError:
            method = ""
        method = method or "1"
    if method == "2":
        nep89 = os.path.join(hq_home(), "examples", "nep", "nep89_20250409.txt")
        if os.path.isfile(nep89):
            shutil.copyfile(nep89, os.path.join(d, "nep.txt"))
            print("nep.txt <- nep89_20250409.txt")
        pc = sorted(glob.glob(os.path.join(d, "POSCAR_*.vasp"))); pc = pc or sorted(glob.glob(os.path.join(d, "*.vasp")))
        if pc:
            import hq_batch
            lat, sp, pos = hq_batch.read_poscar(pc[0])
            hq_batch.write_model_xyz(os.path.join(d, "model.xyz"), lat, sp, pos)
            print(f"model.xyz <- {os.path.basename(pc[0])}")
        with open(os.path.join(d, "run.in"), "w") as f:
            f.write("potential     nep.txt\ncompute_elastic 0.005\n")
        print("run.in <- compute_elastic 0.005 (NEP89)")
    elif method == "3":
        eam_dir = os.path.join(hq_home(), "examples", "elastic", "eam")
        for f in os.listdir(eam_dir):
            s = os.path.join(eam_dir, f)
            if os.path.isfile(s) and not os.path.isfile(os.path.join(d, f)):
                shutil.copyfile(s, os.path.join(d, f))
        print(f"LAMMPS EAM files copied from {eam_dir}")
        print("run: python3 run_wsl_lmp_10queue.py   (or: lmp -in lammps_elastic_eam.in)")
    else:
        print("method 1 (GPUMD 0K): run.in/nep.txt/model.xyz ready")
        print("next: HQmdstemkit.sh 4 402 DIR  (compute_elastic 0K)")
    print(f"elastic inputs ready in {d}: model.xyz, nep.txt, run.in")

def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "strain": cmd_strain(args)
    elif cmd == "0k": cmd_0k(args)
    elif cmd == "plot": cmd_plot(args)
    elif cmd == "born": cmd_born(args)
    elif cmd == "auto": cmd_auto(args)
    elif cmd == "born-plot": cmd_born_plot(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()

