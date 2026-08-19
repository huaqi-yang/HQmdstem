#!/bin/bash
# =============================================================================
# HQmdstemkit Installation  (GPUMDkit style: no big env copy, open-box)
#   - reuses an existing python (nepkit env > HQmdstemkit env > system python3)
#   - no conda env cloning, no massive downloads
#   - writes config.json, adds PATH, links HQmdstemkit.sh into conda base bin
# =============================================================================

echo "======================================================"
echo "  HQmdstemkit Installation"
echo "======================================================"

INSTALL_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo " [1/6] Installing directory: ${INSTALL_DIR}"

# ---- shell rc file ----
RC_FILE="$HOME/.bashrc"
if [[ "$SHELL" == *"zsh"* ]] && [ -f "$HOME/.zshrc" ]; then
    RC_FILE="$HOME/.zshrc"
fi
echo " [2/6] Shell config: ${RC_FILE}"

# ---- locate conda (PATH + common install paths) ----
find_conda() {
    command -v conda 2>/dev/null && return
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
             /opt/miniconda3/bin/conda \
             /opt/anaconda3/bin/conda; do
        [ -x "$p" ] && { echo "$p"; return; }
    done
}
CONDA_BIN="$(find_conda)"
CONDA_ROOT=""
if [ -n "$CONDA_BIN" ]; then
    CONDA_ROOT="$("$CONDA_BIN" info --base 2>/dev/null)"
fi

# ---- choose python with deps: nepkit > HQmdstemkit env > system python3 ----
ENV_PY=""
for c in "${CONDA_ROOT}/envs/nepkit/bin/python3" \
         "${CONDA_ROOT}/envs/HQmdstemkit/bin/python3" \
         "$(command -v python3 2>/dev/null)"; do
    if [ -n "$c" ] && [ -x "$c" ]; then
        if "$c" -c "import numpy, matplotlib, scipy, PIL, pandas" >/dev/null 2>&1; then
            ENV_PY="$c"; break
        fi
        [ -z "$ENV_PY" ] && ENV_PY="$c"
    fi
done
if [ -z "$ENV_PY" ]; then
    echo "[ERROR] python3 not found; install Miniconda first"
    exit 1
fi
echo " [3/6] Using python: ${ENV_PY}"

# ---- ensure python dependencies (only if missing) ----
if "$ENV_PY" -c "import numpy, matplotlib, scipy, PIL, pandas" >/dev/null 2>&1; then
    echo " [4/6] Python dependencies already present (skip install)"
else
    echo " [4/6] Installing python dependencies (numpy matplotlib scipy pillow pandas ase) ..."
    "$ENV_PY" -m ensurepip --upgrade >/dev/null 2>&1 || true
    "$ENV_PY" -m pip install numpy matplotlib scipy pillow pandas ase
fi

# ---- detect external executables ----
find_exe() {
  local name="$1"; shift
  local p
  p="$(command -v "$name" 2>/dev/null)"
  [ -n "$p" ] && { echo "$p"; return; }
  for p in "$@"; do
    [ -x "$p" ] && { echo "$p"; return; }
  done
  echo "$name"
}
GPUMD="$(find_exe gpumd)"
QSTEM="$(find_exe qstem)"
GNEP="$(find_exe gnep)"
NEP="$(find_exe nep)"

# ---- write config.json ----
echo " [5/6] Writing config.json ..."
"$ENV_PY" - "$INSTALL_DIR" "$CONDA_BIN" "$ENV_PY" "$GPUMD" "$QSTEM" "$GNEP" "$NEP" <<'PY'
import json, sys
home, conda_bin, env_py, gpumd, qstem, gnep, nep = sys.argv[1:]
cfg = {
    "env": "nepkit" if "/nepkit/" in env_py else "HQmdstemkit",
    "conda": conda_bin,
    "python": env_py,
    "gpumd": gpumd,
    "qstem": qstem,
    "gnep": gnep,
    "nep": nep,
    "gpumd_example_dir": "examples/gpumd",
    "nep_example_dir": "examples/nep",
    "elastic_example_dir": "examples/elastic",
    "eam_input_dir": "examples/elastic/eam",
    "potential_cuzn_specific": "examples/gpumd/nep.txt",
    "potential_cuzn_generic": "examples/gpumd/nep_gen159000.txt",
    "potential_nep89_universal": "examples/nep/nep89_20250409.txt",
    "potential_nep89_restart": "examples/nep/nep89_20250409.restart",
    "potential_nep89_nepin": "examples/nep/nep89_20250409.nep.in",
    "potential_eam_cuzn": "examples/elastic/eam/CuZn.eam.alloy",
    "train_xyz_example": "examples/nep/train.xyz",
    "select_xyz_example": "examples/nep/selectsum23n16.xyz",
    "rdf_base": "examples/rdf_data",
    "chain_analysis_script": "scripts/hq_tem_core.py",
    "qstem_example_qsc": "examples/qstem/46qstem.qsc",
    "elastic_raw_csv": "examples/elastic_constants_raw_data.csv",
}
json.dump(cfg, open(home + "/config.json", "w"), indent=2, ensure_ascii=False)
print("       config.json written")
PY

# ---- PATH setup (GPUMDkit style) ----
if grep -q "export HQmdstemkit_path=" "$RC_FILE" 2>/dev/null; then
    echo " [6/6] HQmdstemkit config already in ${RC_FILE}"
else
    {
        echo ""
        echo "########### HQmdstemkit Configuration ###########"
        echo "export HQmdstemkit_path=${INSTALL_DIR}"
        echo "export PATH=\${HQmdstemkit_path}:\${PATH}"
        echo "################################################"
    } >> "$RC_FILE"
    echo " [6/6] Added HQmdstemkit PATH to ${RC_FILE}"
fi
chmod +x "${INSTALL_DIR}/HQmdstemkit.sh"

# extra: link into conda base bin so the CURRENT shell can use it right away
if [ -n "$CONDA_ROOT" ]; then
    mkdir -p "$CONDA_ROOT/bin"
    ln -sf "${INSTALL_DIR}/HQmdstemkit.sh" "$CONDA_ROOT/bin/HQmdstemkit.sh"
    chmod +x "$CONDA_ROOT/bin/HQmdstemkit.sh"
    echo "       linked into conda base bin: ${CONDA_ROOT}/bin/HQmdstemkit.sh"
fi

# make variables available in this shell
source "$RC_FILE" 2>/dev/null || true

# ---- self check ----
"$ENV_PY" "$INSTALL_DIR/scripts/hq_env.py" check || true

echo ""
echo "======================================================"
echo "  Installation Complete!  HQmdstemkit is ready to use."
echo "  Usage:  HQmdstemkit.sh"
echo "======================================================"