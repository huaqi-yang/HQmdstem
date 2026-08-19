#!/usr/bin/env python3
"""HQmdstemkit: TEM image analysis (white atomic chain length distribution).

Usage:
  hq_tem.py chain IMAGE [--threshold 150] [--bin-width 0.1] [--output-dir OUT]
"""
import os, sys, shutil, subprocess

USAGE = __doc__

LOCAL_CORE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hq_tem_core.py")
DEFAULT_CORE = [LOCAL_CORE]

def find_core():
    cfg = load_config()
    cp = cfg.get("chain_analysis_script", "")
    if cp and os.path.isfile(cp):
        return cp
    for p in DEFAULT_CORE:
        if os.path.isfile(p):
            return p
    if os.path.isfile(LOCAL_CORE):
        return LOCAL_CORE
    return None

def main():
    if len(sys.argv) < 3 or sys.argv[1] != "chain":
        print(USAGE)
        sys.exit(1)
    image = sys.argv[2]
    if not os.path.isfile(image):
        print(f"image not found: {image}")
        print(USAGE)
        sys.exit(1)
    core = find_core()
    if core is None:
        print("chain_length_analysis.py not found; run from the bondangelGr figures folder")
        print(USAGE)
        sys.exit(1)
    cmd = [sys.executable, core, image] + sys.argv[3:]
    subprocess.run(cmd, check=False)

if __name__ == "__main__":
    main()