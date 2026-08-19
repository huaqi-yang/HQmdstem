#!/usr/bin/env python3
"""Common helpers for HQmdstemkit scripts."""
import json, os, re, sys

def fail(usage):
    print("ERROR: wrong arguments or missing required files.", file=sys.stderr)
    print(usage)
    sys.exit(1)

def require_file(path, usage):
    if not os.path.isfile(path):
        fail(f"required file not found: {path}\n\n{usage}")

def read_extxyz_frames(path):
    """Read extended-xyz frames. Returns list of [n, header_line, atom_lines]."""
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        lines = f.readlines()
    frames = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line.split()[0])
        except ValueError:
            i += 1
            continue
        if i + 1 >= len(lines):
            break
        frames.append([n, lines[i + 1], lines[i + 2:i + 2 + n]])
        i += 2 + n
    return frames

def write_extxyz_frames(frames, path):
    with open(path, "w", encoding="utf-8") as f:
        for n, header, atoms in frames:
            f.write(str(n) + "\n")
            f.write(header if header.endswith("\n") else header + "\n")
            for al in atoms:
                f.write(al if al.endswith("\n") else al + "\n")

def header_energy(header):
    m = re.search(r"energy=([-+0-9.eE]+)", header)
    return float(m.group(1)) if m else None

def set_header_energy(header, energy):
    if re.search(r"energy=([-+0-9.eE]+)", header):
        return re.sub(r"energy=([-+0-9.eE]+)", f"energy={energy:.12g}", header)
    return header.rstrip("\n") + f" energy={energy:.12g}\n"

def frame_key(atoms):
    """Sorted (species, x, y, z) rounded key, insensitive to atom order."""
    pos = []
    for al in atoms:
        c = al.split()
        if len(c) >= 4:
            pos.append((c[0], round(float(c[1]), 2), round(float(c[2]), 2), round(float(c[3]), 2)))
    pos.sort()
    return str(pos)

def parse_thermo_out(path, skip_frames=0):
    """Parse GPUMD thermo.out 18-column rows (T K U P.. box a/b/c vectors)."""
    headers = ["T", "K", "U", "Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy",
               "ax", "ay", "az", "bx", "by", "bz", "cx", "cy", "cz"]
    data = {h: [] for h in headers}
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            parts = line.split()
            if len(parts) != len(headers):
                continue
            try:
                vals = [float(v) for v in parts]
            except ValueError:
                continue
            for h, v in zip(headers, vals):
                data[h].append(v)
    for h in headers:
        data[h] = data[h][skip_frames:]
    return data

def hq_home():
    return os.environ.get("HQMDSTEMKIT_HOME",
                          os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RELATIVE_KEYS = {
    "gpumd_example_dir", "nep_example_dir", "elastic_example_dir", "eam_input_dir",
    "potential_cuzn_specific", "potential_cuzn_generic", "potential_nep89_universal",
    "potential_nep89_restart", "potential_nep89_nepin", "potential_eam_cuzn",
    "train_xyz_example", "select_xyz_example", "abacus_pp_dir", "abacus_orbital_dir",
    "rdf_base", "chain_analysis_script", "qstem_example_qsc", "elastic_raw_csv",
}

def load_config():
    """Load config.json (machine-specific paths); returns dict (possibly empty).
    Package-relative paths in RELATIVE_KEYS are resolved against HQMDSTEMKIT_HOME."""
    cfg = {}
    p = os.path.join(hq_home(), "config.json")
    if os.path.isfile(p):
        try:
            with open(p, encoding="utf-8-sig", errors="replace") as f:
                cfg = json.load(f)
        except Exception:
            cfg = {}
    for k in list(cfg):
        v = cfg.get(k)
        if k in RELATIVE_KEYS and v and not os.path.isabs(v):
            cfg[k] = os.path.join(hq_home(), v)
    return cfg