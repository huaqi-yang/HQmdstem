#!/usr/bin/env python3
"""HQmdstemkit: QSTEM high-resolution TEM workflow (xyz -> cfg -> qstem -> images).

Reference example:
  examples/qstem/
    46qstem.qsc, POSCAR46zhengjiaolammps.cfg, QSTEMdanbaoconvert1.py

Usage:
  hq_qstem.py xyz2cfg IN.xyz [OUT.cfg] [--type Cu=58,Zn=30]
  hq_qstem.py prepare MODEL.xyz OUTDIR [--qsc template.qsc]
  hq_qstem.py convert CFG_FILE [--outdir OUTDIR]
  hq_qstem.py stem IN1.xyz IN2.xyz ... --qsc template.qsc [--outdir OUT] [--ncpu 8] [--scan x0,x1,y0,y1,nx,ny]
  hq_qstem.py run [DIR] [--qstem /path/to/qstem] [--ncpu 8]
  hq_qstem.py list [DIR]
"""
import os, sys, glob, re, subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hq_common import read_extxyz_frames, load_config, hq_home

USAGE = __doc__
_EXAMPLE_QSC_DEFAULT = os.path.join(hq_home(), "examples", "qstem", "46qstem.qsc")

def _cfg():
    return load_config()

def EXAMPLE_QSC():
    return _cfg().get("qstem_example_qsc", _EXAMPLE_QSC_DEFAULT)

def QSTEM_EXE():
    return _cfg().get("qstem", "qstem")
EXAMPLE_CONVERT = os.path.join(hq_home(), "examples", "qstem", "QSTEMdanbaoconvert1.py")
DEFAULT_TYPES = {"Cu": 58, "Zn": 30}

def _copy(src, dst):
    if os.path.isfile(src):
        with open(src, "rb") as a, open(dst, "wb") as b:
            b.write(a.read())
        return True
    return False

def _lattice_from_header(header):
    m = re.search(r'Lattice="([^"]+)"', header)
    if not m:
        return None
    v = [float(x) for x in m.group(1).split()]
    if len(v) < 9:
        return None
    return np.array([v[0:3], v[3:6], v[6:9]])

def _atom_species_pos(atoms):
    sp, pos = [], []
    for al in atoms:
        c = al.split()
        if len(c) >= 4:
            sp.append(c[0])
            pos.append([float(c[1]), float(c[2]), float(c[3])])
    return sp, np.array(pos)

def xyz_to_cfg(xyz_path, cfg_path, type_map=None):
    import numpy as np
    type_map = type_map or DEFAULT_TYPES
    frames = read_extxyz_frames(xyz_path)
    if not frames:
        raise RuntimeError(f"no frame in {xyz_path}")
    n, header, atoms = frames[0]
    lat = _lattice_from_header(header)
    if lat is None:
        raise RuntimeError(f"cannot parse Lattice from {xyz_path}")
    sp, pos = _atom_species_pos(atoms)
    frac = pos @ np.linalg.inv(lat)
    species = sorted(set(sp), key=lambda s: type_map.get(s, 999))
    with open(cfg_path, "w") as f:
        f.write(f"Number of particles = {len(atoms)}\n")
        f.write("A = 1.0 Angstrom (basic length-scale)\n")
        for i in range(3):
            for j in range(3):
                f.write(f"H0({i+1},{j+1}) = {lat[i,j]:.6f} A\n")
        f.write(".NO_VELOCITY.\n")
        f.write(f"entry_count = {len(species)}\n")
        for s in species:
            tid = type_map.get(s, 99)
            f.write(f"{tid}\n{s}\n")
            for k, s0 in enumerate(sp):
                if s0 == s:
                    f.write(f"{frac[k,0]:.8f} {frac[k,1]:.8f} {frac[k,2]:.8f}\n")
    print(f"xyz -> cfg: {cfg_path}  ({len(atoms)} atoms, species {species})")

def _write_qsc_from_template(template, cfg_path, out_qsc, scan=None):
    with open(template, encoding="utf-8-sig", errors="replace") as f:
        text = f.read()
    text = re.sub(r'(?m)^filename:.*$', f'filename: "{cfg_path}"', text)
    if scan:
        x0, x1, y0, y1, nx, ny = scan
        text = re.sub(r'(?m)^scan_x_start:.*$', f"scan_x_start:  {x0:.4f}", text)
        text = re.sub(r'(?m)^scan_x_stop:.*$', f"scan_x_stop:  {x1:.4f}", text)
        text = re.sub(r'(?m)^scan_y_start:.*$', f"scan_y_start:  {y0:.4f}", text)
        text = re.sub(r'(?m)^scan_y_stop:.*$', f"scan_y_stop:  {y1:.4f}", text)
        text = re.sub(r'(?m)^scan_x_pixels:.*$', f"scan_x_pixels: {nx}", text)
        text = re.sub(r'(?m)^scan_y_pixels:.*$', f"scan_y_pixels: {ny}", text)
    with open(out_qsc, "w", encoding="utf-8") as f:
        f.write(text)

def _run_one(exe, qsc):
    d = os.path.dirname(os.path.abspath(qsc)) or "."
    try:
        subprocess.run([exe, os.path.basename(qsc)], cwd=d, check=True, timeout=7200)
        return (qsc, True, "")
    except Exception as e:
        return (qsc, False, str(e))

def _run_all(exe, qscs, ncpu):
    print(f"QSTEM: {len(qscs)} jobs, {ncpu} workers, exe={exe}")
    ok = 0
    with ThreadPoolExecutor(max_workers=ncpu) as pool:
        for qsc, success, err in pool.map(lambda q: _run_one(exe, q), qscs):
            if success:
                ok += 1
                print(f"[OK] {qsc}")
            else:
                print(f"[FAIL] {qsc}: {err}")
    print(f"done: {ok}/{len(qscs)} succeeded")

def cmd_xyz2cfg(args):
    if len(args) < 1:
        print(USAGE); sys.exit(1)
    xyz = args[0]
    cfg = args[1] if len(args) > 1 else os.path.splitext(xyz)[0] + ".cfg"
    type_map = DEFAULT_TYPES
    if "--type" in args:
        i = args.index("--type")
        type_map = {}
        for kv in args[i+1].split(","):
            if "=" in kv:
                s, t = kv.split("=")
                type_map[s] = int(t)
    try:
        xyz_to_cfg(xyz, cfg, type_map)
    except Exception as e:
        print(f"ERROR: {e}\n{USAGE}"); sys.exit(1)

def cmd_prepare(args):
    if len(args) < 2:
        print(USAGE); sys.exit(1)
    model, outdir = args[0], args[1]
    qsc = EXAMPLE_QSC()
    if "--qsc" in args:
        qsc = args[args.index("--qsc") + 1]
    if not os.path.isfile(model):
        print(f"model.xyz not found: {model}\n{USAGE}"); sys.exit(1)
    os.makedirs(outdir, exist_ok=True)
    cfg = os.path.join(outdir, "model.cfg")
    xyz_to_cfg(model, cfg)
    if _copy(qsc, os.path.join(outdir, "model.qsc")):
        _write_qsc_from_template(qsc, cfg, os.path.join(outdir, "model.qsc"))
        print(f"copied/updated QSTEM input -> {os.path.join(outdir, 'model.qsc')}")
    else:
        print("[warn] no QSTEM template found; use --qsc template.qsc")
    _copy(EXAMPLE_CONVERT, os.path.join(outdir, "QSTEMdanbaoconvert1.py"))
    print(f"done. next: HQmdstemkit.sh qstem run {outdir} --ncpu 8")

def cmd_convert(args):
    if len(args) < 1:
        print(USAGE); sys.exit(1)
    cfg = args[0]
    outdir = os.path.dirname(os.path.abspath(cfg)) or "."
    if "--outdir" in args:
        outdir = args[args.index("--outdir") + 1]
    if not os.path.isfile(cfg):
        print(f"cfg file not found: {cfg}\n{USAGE}"); sys.exit(1)
    conv = os.path.join(outdir, "QSTEMdanbaoconvert1.py")
    if os.path.isfile(conv):
        subprocess.run([sys.executable, conv], cwd=outdir, check=False)
        return
    print("QSTEMdanbaoconvert1.py not found; convert manually")

def cmd_stem(args):
    """Full pipeline: xyz files -> cfg -> qsc -> multi-core qstem images."""
    if len(args) < 1:
        print(USAGE); sys.exit(1)
    xyzs = []
    qsc = None
    outdir = "qstem_out"
    ncpu = 8
    exe = QSTEM_EXE()
    scan = None
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--qsc" and i+1 < len(args): qsc = args[i+1]; i += 2
        elif a == "--outdir" and i+1 < len(args): outdir = args[i+1]; i += 2
        elif a == "--ncpu" and i+1 < len(args): ncpu = int(args[i+1]); i += 2
        elif a == "--qstem" and i+1 < len(args): exe = args[i+1]; i += 2
        elif a == "--scan" and i+1 < len(args):
            vals = [float(v) for v in args[i+1].split(",")]
            scan = tuple(int(v) if k >= 4 else v for k, v in enumerate(vals))
            i += 2
        else:
            xyzs.append(a); i += 1
    if not xyzs:
        print("no xyz files given\n" + USAGE); sys.exit(1)
    if not qsc or not os.path.isfile(qsc):
        print(f"QSTEM template qsc not found: {qsc}\n{USAGE}"); sys.exit(1)
    os.makedirs(outdir, exist_ok=True)
    qscs = []
    for xyz in xyzs:
        if not os.path.isfile(xyz):
            print(f"skip missing: {xyz}"); continue
        stem = os.path.splitext(os.path.basename(xyz))[0]
        cfg = os.path.join(outdir, stem + ".cfg")
        q = os.path.join(outdir, stem + ".qsc")
        try:
            xyz_to_cfg(xyz, cfg)
            _write_qsc_from_template(qsc, cfg, q, scan)
            qscs.append(q)
        except Exception as e:
            print(f"[FAIL] {xyz}: {e}")
    _run_all(exe, qscs, ncpu)

def cmd_run(args):
    d = "."
    exe = QSTEM_EXE()
    ncpu = 8
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--qstem" and i+1 < len(args): exe = args[i+1]; i += 2
        elif a == "--ncpu" and i+1 < len(args): ncpu = int(args[i+1]); i += 2
        elif a == "--qstem" or a == "--ncpu": print(USAGE); sys.exit(1)
        else: d = a; i += 1
    qscs = sorted(glob.glob(os.path.join(d, "*.qsc")))
    if not qscs:
        print(f"no *.qsc files in {d}\n{USAGE}"); sys.exit(1)
    _run_all(exe, qscs, ncpu)

def cmd_list(args):
    d = args[0] if args else "."
    for pat in ("*.qsc", "*.cfg", "*.CFG", "QSTEMdanbaoconvert*.py", "*.tif", "*.png"):
        for f in sorted(glob.glob(os.path.join(d, pat))):
            print(f)

def main():
    if len(sys.argv) < 2:
        print(USAGE); sys.exit(1)
    cmd = sys.argv[1]
    args = sys.argv[2:]
    if cmd == "xyz2cfg": cmd_xyz2cfg(args)
    elif cmd == "prepare": cmd_prepare(args)
    elif cmd == "convert": cmd_convert(args)
    elif cmd == "stem": cmd_stem(args)
    elif cmd == "run": cmd_run(args)
    elif cmd == "list": cmd_list(args)
    else: print(USAGE); sys.exit(1)

if __name__ == "__main__":
    main()