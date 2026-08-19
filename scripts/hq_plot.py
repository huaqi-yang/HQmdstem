#!/usr/bin/env python3
"""HQmdstemkit: NEP plotting wrappers (portable, run in any working dir).

Usage:
  hq_plot.py predict DIR [--out nep_prediction_2x2.png]
  hq_plot.py descriptor DIR
  hq_plot.py umap DIR nep.txt ELEMENT POOL.xyz [N_FPS] [train.xyz ...] [--already xyz] [--mode 1|2]
"""
import os, sys, subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, hq_home

USAGE = __doc__
PLOT_DIR = os.path.join(hq_home(), "examples", "plot_scripts")


def _run(script, cwd, args=None):
    p = os.path.join(PLOT_DIR, script)
    if not os.path.isfile(p):
        fail(f"script not found: {p}\n\n{USAGE}")
    os.makedirs(cwd, exist_ok=True)
    subprocess.run([sys.executable, p] + (args or []), cwd=cwd, check=False)


def cmd_predict(args):
    d = args[0] if args else "."
    _run("nep_prediction_beautified.py", d)


def cmd_descriptor(args):
    d = args[0] if args else "."
    _run("nep_descriptor_plots.py", d)


def cmd_umap(args):
    if len(args) < 3:
        fail(USAGE)
    d, nep, elem, pool = args[0], args[1], args[2], args[3]
    _run("umap_select.py", d, [nep, elem, pool] + args[4:])


def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "predict": cmd_predict(args)
    elif cmd == "descriptor": cmd_descriptor(args)
    elif cmd == "umap": cmd_umap(args)
    else: fail(USAGE)


if __name__ == "__main__":
    main()