#!/usr/bin/env python3
"""HQmdstemkit: environment self-check and run results summary.

Usage:
  hq_env.py check [--json]
  hq_env.py config
  hq_env.py summary [DIR]
"""
import glob, importlib, json, os, shutil, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import load_config, hq_home, parse_thermo_out

USAGE = __doc__
REQUIRED_PKGS = ["numpy", "matplotlib", "scipy", "PIL", "pandas"]

def pkg_version(name):
    try:
        m = importlib.import_module(name)
        return getattr(m, "__version__", "?")
    except Exception:
        return None

def exe_ok(path):
    if not path:
        return False
    if str(path).lower().endswith((".exe", ".bat", ".cmd")):
        return os.path.isfile(path)
    if not path:
        return False
    if os.path.sep in path or (os.path.altsep and os.path.altsep in path):
        return os.path.isfile(path) and os.access(path, os.X_OK)
    return shutil.which(path) is not None

def cmd_check(args):
    res = {}
    print("== HQmdstemkit environment check ==")
    print("python:", sys.version.split()[0], sys.executable)
    for pkg in REQUIRED_PKGS:
        v = pkg_version(pkg)
        st = "OK" if v else "MISSING"
        print(f"  {pkg:<12} {v or 'missing':<12} [{st}]")
        res[pkg] = st
    cfg = load_config()
    for key in ("conda", "gpumd", "qstem", "gnep", "nep"):
        exe = cfg.get(key, key)
        ok = exe_ok(exe)
        print(f"  {key:<12} {exe:<24} [{'OK' if ok else 'MISSING'}]")
        res[key] = "OK" if ok else "MISSING"
    print("config:", os.path.join(hq_home(), "config.json"))
    bad = [k for k, v in res.items() if v == "MISSING"]
    print("MISSING:", bad if bad else "none")
    if "--json" in args:
        print(json.dumps(res, indent=2, ensure_ascii=False))

def cmd_config(args):
    cfg = load_config()
    print(json.dumps(cfg, indent=2, ensure_ascii=False))
    print("config file:", os.path.join(hq_home(), "config.json"))

def cmd_summary(args):
    base = args[0] if args else "."
    print(f"== results summary under {base} ==")
    counts = {"model.xyz": 0, "nep.txt": 0, "run.in": 0, "thermo.out": 0,
              "elastic.out": 0, "train.xyz": 0, "nep.in": 0, "movie.xyz": 0,
              "chain_length_stats.csv": 0, "*.qsc": 0, "restart*": 0}
    for root, _, files in os.walk(base):
        for fn in files:
            full = os.path.join(root, fn)
            if fn in ("model.xyz", "nep.txt", "run.in", "thermo.out", "elastic.out",
                      "train.xyz", "nep.in", "movie.xyz", "chain_length_stats.csv"):
                counts[fn] += 1
            elif fn.endswith(".qsc"):
                counts["*.qsc"] += 1
            elif fn.startswith("nep") and fn.endswith(".restart"):
                counts["restart*"] += 1
    for k, v in counts.items():
        print(f"  {k:<22} {v}")
    # key outputs
    for root, _, files in os.walk(base):
        if "elastic.out" in files:
            p = os.path.join(root, "elastic.out")
            try:
                mat = np.loadtxt(p, comments="#")
                print(f"  elastic: {p}  C11={mat[0,0]:.2f} GPa")
            except Exception:
                print(f"  elastic: {p}  (unreadable)")
        if "thermo.out" in files:
            p = os.path.join(root, "thermo.out")
            try:
                d = parse_thermo_out(p, skip_frames=0)
                V = np.mean(np.array(d["ax"])*np.array(d["by"])*np.array(d["cz"]))
                print(f"  thermo : {p}  T~{np.mean(d['T']):.0f} K  V~{V:.1f} A^3")
            except Exception:
                pass
        if "chain_length_stats.csv" in files:
            p = os.path.join(root, "chain_length_stats.csv")
            try:
                import csv
                n = sum(1 for _ in csv.DictReader(open(p, encoding="utf-8-sig")))
                print(f"  chains : {p}  {n} chains")
            except Exception:
                pass

def main():
    if len(sys.argv) < 2:
        print(USAGE); sys.exit(1)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "check": cmd_check(args)
    elif cmd == "config": cmd_config(args)
    elif cmd == "summary": cmd_summary(args)
    else: print(USAGE); sys.exit(1)

if __name__ == "__main__":
    main()