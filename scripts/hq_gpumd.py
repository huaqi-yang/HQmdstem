#!/usr/bin/env python3
"""HQmdstemkit: GPUMD MD simulation (heating / phase transition) workflow.

Example inputs (only model.xyz + run.in + nep.txt are needed):
  examples/gpumd/  (or your own model.xyz + run.in + nep.txt directory)

Usage:
  hq_gpumd.py sample DIR [--mode mttk|scr|mcmd] [--elements Cu,Zn]
  hq_gpumd.py prepare DIR [--T T0] [--dT dT] [--Tmax Tmax] [--steps N] [--dt fs]
  hq_gpumd.py run DIR [--exe gpumd]
  hq_gpumd.py thermo DIR [--out thermo.png]
  hq_gpumd.py phase DIR [--out phase.png]
"""
import os, re, sys, subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, parse_thermo_out, load_config, hq_home

USAGE = __doc__

def ask_choice(prompt, default):
    if not sys.stdin.isatty():
        return default
    try:
        v = input(f"{prompt} [{default}]: ").strip()
    except EOFError:
        v = ""
    return v or default

RUNIN_TEMPLATE = """\
potential   ./nep.txt
velocity    {T}

ensemble    npt_scr {T} {T} 100 0 0 0 0 0 0 100 100 100 100 100 100 50
time_step   {dt}
dump_thermo 100
run         {steps}
"""

MTK_RUNIN = """\
potential   ./nep.txt
velocity    100

ensemble    npt_mttk temp 100 400 aniso 0 0
run         100000

ensemble    npt_mttk temp 400 800 aniso 0 0
dump_thermo 10
dump_exyz   10000
run         100000
"""

SCR_RUNIN = """\
potential   ./nep.txt
velocity    300

ensemble    npt_scr 300 300 100 0 0 0 20 20 100 1000
time_step   1
dump_thermo 1000
dump_position 3000
run         3000000
"""

def cmd_sample(args):
    """GPUMD sampling: MTK heating / NPT_SCR / MC-MD (choose elements & potential)."""
    d = None; mode = None; elements = None
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--mode" and i+1 < len(args): mode = args[i+1].lower(); i += 2
        elif a == "--elements" and i+1 < len(args): elements = args[i+1]; i += 2
        elif not a.startswith("-") and d is None: d = a; i += 1
        else: fail(f"unknown option {a}\n\n{USAGE}")
    d = d or os.getcwd()
    os.makedirs(d, exist_ok=True)
    if mode is None:
        mode = ask_choice("GPUMD sampling mode (mttk/scr/mcmd)", "mttk").lower()
    if mode not in ("mttk", "scr", "mcmd"):
        fail(f"--mode must be mttk, scr or mcmd\n\n{USAGE}")
    if elements is None:
        elements = ask_choice("Element species (comma separated, e.g. Cu,Zn)", "Cu,Zn")
    sel = set(e.capitalize() for e in re.split(r"[,，\s]+", elements) if e)
    is_cuzn = {"Cu", "Zn"} <= sel
    if not os.path.isfile(os.path.join(d, "model.xyz")):
        src = load_config().get("gpumd_example_dir") or os.path.join(hq_home(), "examples", "gpumd")
        with open(os.path.join(src, "model.xyz"), encoding="utf-8-sig", errors="replace") as a, \
             open(os.path.join(d, "model.xyz"), "w", encoding="utf-8") as b:
            b.write(a.read())
    if is_cuzn:
        pot_src = os.path.join(hq_home(), "examples", "gpumd", "nep.txt")
        print(f"[GPUMD] elements {','.join(sorted(sel))} -> CuZn dedicated potential (nep.txt)")
    else:
        pot_src = os.path.join(hq_home(), "examples", "nep", "nep89_20250409.txt")
        print(f"[GPUMD] elements {','.join(sorted(sel))} -> generic NEP89 (copied as nep.txt)")
    if os.path.isfile(pot_src):
        with open(pot_src, encoding="utf-8-sig", errors="replace") as a, \
             open(os.path.join(d, "nep.txt"), "w", encoding="utf-8") as b:
            b.write(a.read())
    if mode == "mttk":
        content = MTK_RUNIN
        label = "MTK heating 100K->400K->800K (NPT_MTTK)"
    elif mode == "scr":
        content = SCR_RUNIN
        label = "NPT_SCR isothermal sampling 300 K"
    else:
        mcmd = os.path.join(hq_home(), "examples", "gpumd", "run_mcmd.in")
        if not os.path.isfile(mcmd):
            fail(f"run_mcmd.in not found: {mcmd}\n\n{USAGE}")
        with open(mcmd, encoding="utf-8-sig", errors="replace") as f:
            content = f.read()
        label = "MC/MD hybrid sampling (20-400 K)"
    with open(os.path.join(d, "run.in"), "w", encoding="utf-8") as f:
        f.write(content)
    print(f"[OK] {label}")
    print(f"[OK] {os.path.join(d, 'run.in')}  model.xyz  nep.txt")
    print("next: HQmdstemkit.sh 3 302 DIR   (or: cd DIR && gpumd)")

def cmd_auto(args):
    """Auto: detect local model.xyz+nep.txt, write run.in; else use bundled example."""
    d = None
    T, steps, dt = 300, 1100000, 1
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--T" and i+1 < len(args): T = int(args[i+1]); i += 2
        elif a == "--steps" and i+1 < len(args): steps = int(args[i+1]); i += 2
        elif a == "--dt" and i+1 < len(args): dt = int(args[i+1]); i += 2
        elif not a.startswith("-") and d is None: d = a; i += 1
        else: fail(f"unknown option {a}\n\n{USAGE}")
    d = d or os.getcwd()
    os.makedirs(d, exist_ok=True)
    if not (os.path.isfile(os.path.join(d, "model.xyz")) and os.path.isfile(os.path.join(d, "nep.txt"))):
        src = load_config().get("gpumd_example_dir") or os.path.join(hq_home(), "examples", "gpumd")
        for f in ("model.xyz", "nep.txt"):
            if not os.path.isfile(os.path.join(src, f)):
                fail(f"example file missing: {os.path.join(src, f)}\n\n{USAGE}")
            with open(os.path.join(src, f), encoding="utf-8-sig", errors="replace") as a, \
                 open(os.path.join(d, f), "w", encoding="utf-8") as b:
                b.write(a.read())
        print(f"no local model.xyz/nep.txt -> copied bundled example from {src}")
    run_path = os.path.join(d, "run.in")
    with open(run_path, "w") as f:
        f.write(RUNIN_TEMPLATE.format(T=T, steps=steps, dt=dt))
    print(f"[OK] {run_path} written (T={T} K, steps={steps}, dt={dt} fs)")
    print("next: HQmdstemkit.sh 3 302   (or: cd DIR && gpumd)")

def cmd_example(args):
    d = args[0] if args else os.path.join(hq_home(), "examples", "gpumd")
    os.makedirs(d, exist_ok=True)
    src = load_config().get("gpumd_example_dir") or os.path.join(hq_home(), "examples", "gpumd")
    for f in ("model.xyz", "nep.txt", "run.in"):
        s = os.path.join(src, f)
        if os.path.isfile(s):
            with open(s, encoding="utf-8-sig", errors="replace") as a, \
                 open(os.path.join(d, f), "w", encoding="utf-8") as b:
                b.write(a.read())
    print(f"example GPUMD inputs ready in {d}: model.xyz, nep.txt, run.in")

def cmd_prepare(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    os.makedirs(d, exist_ok=True)
    for f in ("model.xyz", "nep.txt"):
        if not os.path.isfile(os.path.join(d, f)):
            fail(f"missing {f} in {d}\n\n{USAGE}")
    T0 = 300; dT = 100; Tmax = 600; steps = 1100000; dt = 1
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--T" and i+1 < len(args): T0 = int(args[i+1]); i += 2
        elif a == "--dT" and i+1 < len(args): dT = int(args[i+1]); i += 2
        elif a == "--Tmax" and i+1 < len(args): Tmax = int(args[i+1]); i += 2
        elif a == "--steps" and i+1 < len(args): steps = int(args[i+1]); i += 2
        elif a == "--dt" and i+1 < len(args): dt = int(args[i+1]); i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    with open(os.path.join(d, "run.in"), "w") as f:
        f.write(RUNIN_TEMPLATE.format(T=T0, steps=steps, dt=dt))
    # heating driver: one folder per temperature
    heatsh = os.path.join(d, "run_heating.sh")
    with open(heatsh, "w") as f:
        f.write("#!/usr/bin/env bash\nset -e\n")
        for T in range(T0, Tmax + 1, dT):
            sub = os.path.join(d, f"{T}K")
            os.makedirs(sub, exist_ok=True)
            with open(os.path.join(sub, "run.in"), "w") as f2:
                f2.write(RUNIN_TEMPLATE.format(T=T, steps=steps, dt=dt))
            with open(os.path.join(d, "model.xyz"), encoding="utf-8-sig", errors="replace") as src, \
                 open(os.path.join(sub, "model.xyz"), "w", encoding="utf-8") as f2:
                f2.write(src.read())
            with open(os.path.join(d, "nep.txt"), encoding="utf-8-sig", errors="replace") as src, \
                 open(os.path.join(sub, "nep.txt"), "w", encoding="utf-8") as f2:
                f2.write(src.read())
            f.write(f"(cd {sub} && gpumd)\n")
    os.chmod(heatsh, 0o755)
    print(f"prepared run.in (T0={T0}, dT={dT}, Tmax={Tmax}, steps={steps}) in {d}")
    print(f"heating driver -> {heatsh}")

def cmd_run(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    exe = "gpumd"
    if len(args) > 1 and args[1] == "--exe":
        exe = args[2]
    for f in ("model.xyz", "run.in", "nep.txt"):
        if not os.path.isfile(os.path.join(d, f)):
            fail(f"missing {f} in {d}\n\n{USAGE}")
    subprocess.run([exe], cwd=d, check=True)
    print("GPUMD finished")

def cmd_thermo(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    tp = os.path.join(d, "thermo.out")
    if not os.path.isfile(tp):
        fail(f"thermo.out not found in {d}\n\n{USAGE}")
    data = parse_thermo_out(tp, skip_frames=0)
    T = data["T"]
    V = [a*b*c for a, b, c in zip(data["ax"], data["by"], data["cz"])]
    print(f"frames={len(T)}  T_mean={sum(T)/len(T):.1f} K  "
          f"V_mean={sum(V)/len(V):.3f} A^3  "
          f"P_mean={sum(data['Pxx'])/len(data['Pxx']):.3f} GPa")
    out = os.path.join(d, "thermo.png")
    if len(args) > 1 and args[1] == "--out":
        out = args[2]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.3, 3.0), dpi=200)
        ax.plot(T, V, lw=1.2)
        ax.set_xlabel("T (K)"); ax.set_ylabel("V ($\\AA^3$)")
        ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)
        print(f"thermo plot -> {out}")
    except Exception as e:
        print(f"[warn] plotting skipped: {e}")

def cmd_phase(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    rows = []
    for name in sorted(os.listdir(d)):
        if name.endswith("K") and os.path.isfile(os.path.join(d, name, "thermo.out")):
            data = parse_thermo_out(os.path.join(d, name, "thermo.out"), skip_frames=100)
            V = [a*b*c for a, b, c in zip(data["ax"], data["by"], data["cz"])]
            rows.append((int(name[:-1]), sum(V)/len(V)))
    if not rows:
        fail(f"no <T>K/thermo.out found in {d}\n\n{USAGE}")
    rows.sort()
    for T, V in rows:
        print(f"{T} K   V = {V:.3f} A^3")
    out = os.path.join(d, "phase.png")
    if len(args) > 1 and args[1] == "--out":
        out = args[2]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        Ts = [r[0] for r in rows]; Vs = [r[1] for r in rows]
        fig, ax = plt.subplots(figsize=(4.3, 3.0), dpi=200)
        ax.plot(Ts, Vs, "o-", lw=1.4, ms=4)
        ax.set_xlabel("T (K)"); ax.set_ylabel("V ($\\AA^3$)")
        ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout(); fig.savefig(out, dpi=200); plt.close(fig)
        print(f"phase plot -> {out}")
    except Exception as e:
        print(f"[warn] plotting skipped: {e}")

def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "sample": cmd_sample(args)
    elif cmd == "auto": cmd_auto(args)
    elif cmd == "example": cmd_example(args)
    elif cmd == "prepare": cmd_prepare(args)
    elif cmd == "run": cmd_run(args)
    elif cmd == "thermo": cmd_thermo(args)
    elif cmd == "phase": cmd_phase(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()