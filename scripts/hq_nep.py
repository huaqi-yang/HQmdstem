#!/usr/bin/env python3
"""HQmdstemkit: NEP training / fine-tuning (nepkit) and dataset tools.

Examples:
  hq_nep.py train DIR --train-xyz train.xyz [--nep-in nep.in] [--exe gnep]
  hq_nep.py finetune DIR [--restart nep89_20250409.restart] [--type "2 Cu Zn"]
  hq_nep.py shift IN.xyz OUT.xyz [--ref E_Cu,E_Zn,E_C,E_Ti,E_Al,E_V]
  hq_nep.py remove SRC.xyz RM.xyz -o OUT.xyz
"""
import os, sys, subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import fail, read_extxyz_frames, write_extxyz_frames, hq_home, \
    header_energy, set_header_energy, frame_key, load_config

USAGE = __doc__

def _example_dir():
    cfg = load_config()
    return cfg.get("nep_example_dir") or os.path.join(hq_home(), "examples", "nep")

def _ensure_example(name, out_dir, out_name):
    out = os.path.join(out_dir, out_name)
    if os.path.isfile(out):
        return False
    src = os.path.join(_example_dir(), name)
    if not os.path.isfile(src):
        return False
    with open(src, encoding="utf-8-sig", errors="replace") as a, \
         open(out, "w", encoding="utf-8") as b:
        b.write(a.read())
    print(f"copied example {name} -> {out}")
    return True

DEFAULT_REFS = {
    "Cu": -4965.389826, "Zn": -5444.916537, "C": -155.0907296,
    "Ti": -1562.941684, "Al": -59.386039, "V": -2000.204008,
}

def cmd_train(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    os.makedirs(d, exist_ok=True)
    train_xyz = os.path.join(d, "train.xyz")
    nep_in = os.path.join(d, "nep.in")
    exe = "gnep"
    i = 1
    while i < len(args):
        a = args[i]
        if a == "--train-xyz" and i+1 < len(args): train_xyz = args[i+1]; i += 2
        elif a == "--nep-in" and i+1 < len(args): nep_in = args[i+1]; i += 2
        elif a == "--exe" and i+1 < len(args): exe = args[i+1]; i += 2
        else: fail(f"unknown option {a}\n\n{USAGE}")
    if not os.path.isfile(train_xyz):
        _ensure_example("selectsum23n16.xyz", d, "train.xyz")
    if not os.path.isfile(train_xyz):
        fail(f"train.xyz not found: {train_xyz}\n\n{USAGE}")
    if not os.path.isfile(nep_in):
        # minimal Cu-Zn NEP4-ZBL template (edit as needed)
        with open(nep_in, "w") as f:
            f.write("nep4_zbl 2 Cu Zn\nzbl 1 2\ncutoff 6 5 97 61\n"
                    "n_max 4 4\nbasis_size 8 8\nl_max 4 2 1\nANN 80 0\n")
    try:
        subprocess.run([exe, "nep.in"], cwd=d, check=True)
    except FileNotFoundError:
        fail(f"executable not found: {exe}\nplease activate the nepkit/HQmdstemkit env or install nep\n\n{USAGE}")
    print("NEP training finished")

def _copy_file(src, dst):
    with open(src, "rb") as a, open(dst, "wb") as b:
        b.write(a.read())

def cmd_finetune(args):
    if len(args) < 1:
        fail(USAGE)
    d = args[0]
    restart = None
    type_line = None
    i = 1
    while i < len(args):
        if args[i] == "--restart" and i+1 < len(args):
            restart = args[i+1]; i += 2
        elif args[i] == "--type" and i+1 < len(args):
            type_line = args[i+1]; i += 2
        else:
            fail(f"unknown option {args[i]}\n\n{USAGE}")
    # training set + final potential
    _ensure_example("selectsum23n16.xyz", d, "train.xyz")
    _ensure_example("nep_best3n16v3.txt", d, "nep.txt")
    # bundled NEP89 files
    exdir = _example_dir()
    for name in ("nep89_20250409.nep.in", "nep89_20250409.restart", "nep89_20250409.txt"):
        _ensure_example(name, d, name)
    # nep.in from the NEP89 fine-tune template
    nep_in = os.path.join(d, "nep.in")
    src = os.path.join(exdir, "nep89_20250409.nep.in")
    if os.path.isfile(src):
        with open(src, encoding="utf-8-sig", errors="replace") as a, \
             open(nep_in, "w", encoding="utf-8") as b:
            b.write(a.read())
        print("nep89_20250409.nep.in -> nep.in")
    # type line: ask the user, e.g. "2 Cu Zn"
    if type_line is None:
        try:
            type_line = input("Input element types (e.g., 2 Cu Zn) [2 Cu Zn]: ").strip()
        except EOFError:
            type_line = ""
        if not type_line:
            type_line = "2 Cu Zn"
    # restart handling
    restart_ref = "nep89_20250409.restart"
    if restart and os.path.isfile(restart):
        _copy_file(restart, os.path.join(d, os.path.basename(restart)))
        restart_ref = os.path.basename(restart)
    if not os.path.isfile(os.path.join(d, restart_ref)):
        fail(f"restart file not found: {os.path.join(d, restart_ref)}\n\n{USAGE}")
    # patch nep.in: type + fine_tune line
    with open(nep_in, encoding="utf-8") as f:
        txt = f.read()
    if "type <your types>" in txt:
        txt = txt.replace("type <your types>", f"type {type_line}")
    lines = []
    ft_written = False
    for ln in txt.splitlines():
        if ln.lstrip().startswith("fine_tune "):
            lines.append(f"fine_tune nep89_20250409.txt {restart_ref}")
            ft_written = True
        else:
            lines.append(ln)
    if not ft_written:
        lines.insert(0, f"fine_tune nep89_20250409.txt {restart_ref}")
    with open(nep_in, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"finetune prepared in {d}:")
    print("  train.xyz (selectsum23n16.xyz)")
    print("  nep.txt   (nep_best3n16v3.txt)")
    print(f"  nep.in    (type {type_line})")
    print(f"  {restart_ref}, nep89_20250409.txt")
    print("run: (cd DIR && gnep nep.in)   in the nepkit/HQmdstemkit environment")

def cmd_shift(args):
    if len(args) < 2:
        fail(USAGE)
    src, out = args[0], args[1]
    refs = dict(DEFAULT_REFS)
    if "--ref" in args:
        idx = args.index("--ref")
        vals = args[idx+1].split(",")
        for key, val in zip(("Cu", "Zn", "C", "Ti", "Al", "V"), vals):
            try: refs[key] = float(val)
            except ValueError: pass
    frames = read_extxyz_frames(src)
    n_ok = 0
    for frame in frames:
        n, header, atoms = frame
        e = header_energy(header)
        if e is None:
            continue
        counts = {}
        for al in atoms:
            c = al.split()
            if c: counts[c[0]] = counts.get(c[0], 0) + 1
        shift = sum(refs.get(sp, 0.0) * cnt for sp, cnt in counts.items())
        frame[1] = set_header_energy(header, e - shift)
        n_ok += 1
    write_extxyz_frames(frames, out)
    print(f"shifted {n_ok} frames -> {out}")

def cmd_remove(args):
    if len(args) < 2:
        fail(USAGE)
    src, rm = args[0], args[1]
    out = "selectsum.xyz"
    if "-o" in args:
        out = args[args.index("-o") + 1]
    src_f = read_extxyz_frames(src)
    rm_f = read_extxyz_frames(rm)
    rm_keys = {frame_key(a) for _, _, a in rm_f}
    kept = []
    found = 0
    for n, h, a in src_f:
        if frame_key(a) in rm_keys:
            found += 1
        else:
            kept.append([n, h, a])
    write_extxyz_frames(kept, out)
    print(f"src={len(src_f)}  remove={len(rm_f)}  match={found}  keep={len(kept)} -> {out}")

def main():
    if len(sys.argv) < 2:
        fail(USAGE)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "train": cmd_train(args)
    elif cmd == "finetune": cmd_finetune(args)
    elif cmd == "shift": cmd_shift(args)
    elif cmd == "remove": cmd_remove(args)
    else: fail(USAGE)

if __name__ == "__main__":
    main()