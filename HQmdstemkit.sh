#!/usr/bin/env bash
# =============================================================================
# HQmdstemkit  (Cu-Zn GPUMD / NEP / QSTEM workflow toolkit)
# Main executable: HQmdstemkit.sh     (menu style follows GPUMDkit)
# =============================================================================

set -u

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" >/dev/null 2>&1 && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [ "${SOURCE#/}" != "$SOURCE" ] || SOURCE="$DIR/$SOURCE"
done
HQMDSTEMKIT_HOME="$(cd -P "$(dirname "$SOURCE")" >/dev/null 2>&1 && pwd)"
SCRIPT_DIR="${HQMDSTEMKIT_HOME}/scripts"
PYTHON="${PYTHON:-python3}"
if [ -f "${HQMDSTEMKIT_HOME}/config.json" ]; then
  CFG_PY="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("python",""))' "${HQMDSTEMKIT_HOME}/config.json" 2>/dev/null)"
  [ -n "${CFG_PY}" ] && PYTHON="${CFG_PY}"
fi
if [ -n "${PYTHON:-}" ] && ! [ -x "${PYTHON}" ] && ! command -v "${PYTHON}" >/dev/null 2>&1; then
  PYTHON="python3"
fi

banner() {
cat <<'EOF'
H   H   QQQ   M   M  DDDD    SSS   TTTTT  EEEEE  M   M  K   K  IIIII  TTTTT
H   H  Q   Q  MM MM  D   D  S        T    E      MM MM  K  K     I      T
HHHHH  Q   Q  M M M  D   D   SSS     T    EEE    M M M  KKK      I      T
H   H  Q  QQ  M   M  D   D      S    T    E      M   M  K  K     I      T
H   H   QQQQ  M   M  DDDD   SSSS     T    EEEEE  M   M  K   K  IIIII    T

     HQmdstemkit (Cu-Zn GPUMD / NEP / QSTEM / CP2K / ABACUS)
EOF
}

ask()    { printf "%s: " "$1"; read -r "$2"; }
run_py() { local s="$1"; shift; "${PYTHON}" "${SCRIPT_DIR}/${s}" "$@"; }
pause()  { read -r -p "Press Enter to continue... " _; }

# ---------------- handlers ----------------
run_xyz2cfg() {
  local in out
  in="${1:-}"; out="${2:-}"
  [ -z "$in" ] && ask "input xyz file" in
  [ -z "$out" ] && ask "output cfg file" out
  [ -z "$out" ] && out="${in%.xyz}.cfg"
  run_py hq_qstem.py xyz2cfg "$in" "$out"
}

run_cfg2xyz() {
  local in
  in="${1:-}"
  [ -z "$in" ] && ask "input cfg file" in
  run_py hq_qstem.py convert "$in"
}

run_shift() {
  local in out
  in="${1:-}"; out="${2:-}"
  [ -z "$in" ] && ask "input xyz (train.xyz)" in
  [ -z "$out" ] && ask "output xyz (shifted.xyz)" out
  run_py hq_nep.py shift "$in" "$out"
}

run_remove() {
  local src rm out
  src="${1:-}"; rm="${2:-}"; out="${3:-}"
  [ -z "$src" ] && ask "source xyz" src
  [ -z "$rm" ]  && ask "remove xyz" rm
  [ -z "$out" ] && out="selectsum.xyz"
  run_py hq_nep.py remove "$src" "$rm" -o "$out"
}

run_structure() {
  local kind n
  kind="${1:-}"; n="${2:-}"
  [ -z "$kind" ] && ask "structure type (L12/B2/L10)" kind
  [ -z "$n" ] && ask "supercell size N" n
  run_py hq_structure.py ordered "$kind" "$n"
}

run_disordered() {
  local ncu nzn
  ncu="${1:-}"; nzn="${2:-}"
  [ -z "$ncu" ] && ask "number of Cu atoms" ncu
  [ -z "$nzn" ] && ask "number of Zn atoms" nzn
  run_py hq_structure.py disordered "$ncu" "$nzn"
}

run_gpumd_auto()   { run_py hq_gpumd.py auto "$@"; }
run_gpumd_sample() { run_py hq_gpumd.py sample "$@"; }
run_gpumd_prepare() {
  local d
  d="${1:-}"
  [ -z "$d" ] && ask "GPUMD directory (model.xyz + nep.txt)" d
  run_py hq_gpumd.py prepare "$d" --T 300 --dT 100 --Tmax 600 --steps 1100000
}

run_gpumd_run() {
  local d
  d="${1:-}"
  [ -z "$d" ] && ask "GPUMD directory" d
  run_py hq_gpumd.py run "$d"
}

run_nep_train() {
  local d
  d="${1:-}"
  [ -z "$d" ] && ask "training directory (train.xyz)" d
  run_py hq_nep.py train "$d"
}

run_nep_finetune() {
  local d r
  d="${1:-}"; r="${2:-}"
  [ -z "$d" ] && ask "training directory" d
  [ -z "$r" ] && ask "restart file (nep89_20250409.restart)" r
  run_py hq_nep.py finetune "$d" --restart "$r"
}

run_qstem_stem() {
  local qsc out
  qsc="${1:-}"; out="${2:-}"
  [ -z "$qsc" ] && ask "QSTEM template qsc" qsc
  [ -z "$out" ] && out="qstem_out"
  run_py hq_qstem.py stem . --qsc "$qsc" --outdir "$out" --ncpu 8
}

run_elastic_0k() {
  local d
  d="${1:-}"
  [ -z "$d" ] && ask "directory (model.xyz + nep.txt)" d
  run_py hq_elastic.py 0k "$d" --strain 0.01 --run
}

run_elastic_strain() {
  local d
  d="${1:-}"
  [ -z "$d" ] && ask "directory with thermo.out" d
  run_py hq_elastic.py strain "$d" --T 300 --skip 1000 --slices 10
}

run_elastic_born() {
  local f
  f="${1:-}"
  [ -z "$f" ] && ask "elastic.out file" f
  run_py hq_elastic.py born "$f"
}

run_hull()      { run_py hq_phase.py hull "${1:-}"; }
run_hull_dft()  { run_py hq_phase.py hull-dft "${1:-}"; }
run_phase_exp() { run_py hq_phase.py exp "${1:-}"; }

run_rdf_4x1() {
  local base
  base="${1:-}"
  [ -z "$base" ] && base="${HQMDSTEMKIT_HOME}/examples/rdf_data"
  run_py hq_rdf.py 4x1 --base "$base"
}

run_rdf_plot() {
  local f
  f="${1:-}"
  [ -z "$f" ] && ask "rdf file (rdf.out / gr.txt)" f
  run_py hq_rdf.py rdf "$f"
}

run_born_plot() {
  local csv
  csv="${1:-}"
  [ -z "$csv" ] && csv="${HQMDSTEMKIT_HOME}/examples/elastic_constants_raw_data.csv"
  run_py hq_elastic.py born-plot --csv "$csv"
}

run_elastic_plot() {
  local csv
  csv="${1:-}"
  [ -z "$csv" ] && csv="${HQMDSTEMKIT_HOME}/examples/elastic_constants_raw_data.csv"
  run_py hq_elastic.py plot "$csv"
}

run_micro() {
  local cmd f
  cmd="${1:-}"; f="${2:-}"
  [ -z "$cmd" ] && ask "micro command (chain/cluster/segregate/twin/grain/orientation)" cmd
  [ -z "$f" ] && ask "input xyz / image file" f
  run_py hq_micro.py "$cmd" "$f"
}

run_thermo() { local d; d="${1:-}"; [ -z "$d" ] && ask "GPUMD directory (thermo.out)" d; run_py hq_gpumd.py thermo "$d"; }
run_phase()  { local d; d="${1:-}"; [ -z "$d" ] && ask "directory with T K subdirs" d; run_py hq_gpumd.py phase "$d"; }
run_qstem_prepare() {
  local m o
  m="${1:-}"; o="${2:-}"
  [ -z "$m" ] && ask "model xyz" m
  [ -z "$o" ] && ask "output dir" o
  run_py hq_qstem.py prepare "$m" "$o"
}
run_qstem_run()   { local d; d="${1:-}"; [ -z "$d" ] && d="."; run_py hq_qstem.py run "$d" --ncpu 8; }
run_qstem_list()  { run_py hq_qstem.py list "${1:-.}"; }
run_env_check()   { run_py hq_env.py check; }
run_elastic_auto() { run_py hq_elastic.py auto "$@"; }
run_nep_predict()    { run_py hq_plot.py predict "$@"; }
run_nep_descriptor() { run_py hq_plot.py descriptor "$@"; }
run_umap()           { run_py hq_plot.py umap "$@"; }
run_cohesive()       { run_py hq_calc.py cohesive "$@"; }
run_shear()          { run_py hq_calc.py shear "$@"; }
run_sf()             { run_py hq_calc.py stacking-fault "$@"; }

run_batch()       { if [ $# -eq 0 ]; then run_py hq_batch.py batch "."; else run_py hq_batch.py batch "$@"; fi; }
run_abacus_pretreat() { if [ $# -eq 0 ]; then run_py hq_abacus.py pretreat "."; else run_py hq_abacus.py pretreat "$@"; fi; }
run_abacus_menu() {
  local f dn
  printf "ABACUS functional? 1) LDA  2) PBE  3) PBE-D3 [3]: "
  read -r f
  case "$f" in
    1|LDA|lda) f="LDA";;
    2|PBE|pbe) f="PBE";;
    3|PBE-D3|pbe-d3|"") f="PBE-D3";;
    *) f="PBE-D3";;
  esac
  ask "DFT data folder name under \$HOME [CuZn]" dn
  [ -z "$dn" ] && dn="CuZn"
  run_py hq_abacus.py pretreat "." --func "$f" --data-name "$dn"
}
run_abacus_extract()  { if [ $# -eq 0 ]; then run_py hq_abacus.py extract "."; else run_py hq_abacus.py extract "$@"; fi; }
run_abacus_shift()    { if [ $# -eq 0 ]; then local f; ask "total_train.xyz" f; run_py hq_abacus.py shift "$f"; else run_py hq_abacus.py shift "$@"; fi; }
run_cp2k_pretreat()   { if [ $# -eq 0 ]; then run_py hq_cp2k.py pretreat "."; else run_py hq_cp2k.py pretreat "$@"; fi; }
run_cp2k_extract()    { if [ $# -eq 0 ]; then run_py hq_cp2k.py extract "."; else run_py hq_cp2k.py extract "$@"; fi; }
run_cp2k_menu() {
  local f dn
  printf "CP2K functional? 1) LDA  2) PBE  3) PBE-D3 [3]: "
  read -r f
  case "$f" in
    1|LDA|lda) f="LDA";;
    2|PBE|pbe) f="PBE";;
    3|PBE-D3|pbe-d3|"") f="PBE-D3";;
    *) f="PBE-D3";;
  esac
  ask "DFT data folder name under \$HOME [CuZn]" dn
  [ -z "$dn" ] && dn="CuZn"
  run_py hq_cp2k.py pretreat "." --func "$f" --data-name "$dn"
}
run_summary()     { local d; d="${1:-}"; [ -z "$d" ] && d="."; run_py hq_env.py summary "$d"; }
show_readme()     { sed -n '1,200p' "${HQMDSTEMKIT_HOME}/README.md" 2>/dev/null || cat "${HQMDSTEMKIT_HOME}/README.md"; }

# ---------------- submenus ----------------
box() { echo " +-------------------------------------------------------------+"; }

menu_format() {
  while true; do
    box
    echo " |                   FORMAT CONVERSION TOOLS                   |"
    box
    echo " | 101) xyz to QSTEM cfg            104) Remove bad frames    |"
    echo " | 102) QSTEM cfg to xyz            105) Replicate structure  |"
    echo " | 103) ABACUS energy shift                                   |"
    box
    echo " | xyz2cfg) xyz to cfg   cfg2xyz) cfg to xyz   shift) shift   |"
    echo " | remove)  remove frames                                     |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number or keyword: "
    read -r opt
    case "$opt" in
      101|xyz2cfg) run_xyz2cfg;;
      102|cfg2xyz) run_cfg2xyz;;
      103|shift)   run_shift;;
      104|remove)  run_remove;;
      105|replicate) echo "replicate: not implemented yet; see hq_structure.py";;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_structures() {
  while true; do
    box
    echo " |                   SAMPLE STRUCTURES TOOLS                 |"
    box
    echo " | 201) L12 Cu3Zn               204) Disordered Cu80Zn20     |"
    echo " | 202) B2 CuZn                 205) Disordered custom       |"
    echo " | 203) L10 CuZn                                              |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number: "
    read -r opt
    case "$opt" in
      201|l12) run_structure L12 2;;
      202|b2)  run_structure B2 2;;
      203|l10) run_structure L10 2;;
      204|cu80zn20) run_disordered 80 20;;
      205|disordered) run_disordered;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_workflow() {
  while true; do
    box
    echo " |                      WORKFLOW TOOLS                       |"
    box
    echo " | 300) GPUMD sampling (MTK/SCR/MCMD) 305) QSTEM full pipeline|"
    echo " | 301) GPUMD prepare (auto/example) 306) Elastic 0K compute   |"
    echo " | 302) GPUMD run                                              |"
    echo " | 303) NEP train                 308) Batch pretreatment      |"
    echo " | 304) NEP finetune              309) ABACUS SCF pretreat     |"
    echo " |                             310) ABACUS extract             |"
    echo " |                             311) ABACUS energy shift        |"
    echo " |                             312) CP2K SCF pretreat           |"
    echo " |                             313) CP2K extract energy         |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number: "
    read -r opt
    case "$opt" in
      300|gpumdsample) run_gpumd_sample;;
      301|auto)     run_gpumd_auto;;
      307|prepare)  run_gpumd_prepare;;
      302|run)      run_gpumd_run;;
      303|train)    run_nep_train;;
      304|finetune) run_nep_finetune;;
      305|stem)     run_qstem_stem;;
      306|elastic0k) run_elastic_0k;;
      308|batch)     run_batch;;
      309|abacus)    run_abacus_menu;;
      310|abacusextract) run_abacus_extract;;
      311|abacusshift)   run_abacus_shift;;
      312|cp2k)          run_cp2k_menu;;
      313|cp2kextract)   run_cp2k_extract;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_calculators() {
  while true; do
    box
    echo " |                      CALCULATOR TOOLS                     |"
    box
    echo " | 401) Elastic auto (0K 1/2/3)    404) Convex hull DFT    |"
    echo " | 402) Elastic 0K                    405) Phase diagram      |"
    echo " | 403) Convex hull                   406) Born stability     |"
    echo " | 408) Cohesive/EOS                409) Shear (todo)         |"
    echo " | 410) Stacking fault (todo)                                 |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number: "
    read -r opt
    case "$opt" in
      401|auto)      run_elastic_auto;;
      407|strain)    run_elastic_strain;;
      402|elastic0k) run_elastic_0k;;
      403|hull)      run_hull;;
      404|hull-dft)  run_hull_dft;;
      405|exp)       run_phase_exp;;
      406|born)      run_elastic_born;;
      408|cohesive)  run_cohesive;;
      409|shear)     run_shear;;
      410|stackingfault) run_sf;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_analyzer() {
  while true; do
    box
    echo " |                       ANALYZER TOOLS                      |"
    box
    echo " | 501) TEM chain length        505) Grain size              |"
    echo " | 502) Cluster analysis        506) Orientation             |"
    echo " | 503) Segregation profile     507) GPUMD thermo            |"
    echo " | 504) Twin detection          508) GPUMD phase             |"
    box
    echo " | chain) cluster) segregate) twin) grain) orientation)       |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number or keyword: "
    read -r opt
    case "$opt" in
      501|chain)       run_micro chain;;
      502|cluster)     run_micro cluster;;
      503|segregate)   run_micro segregate;;
      504|twin)        run_micro twin;;
      505|grain)       run_micro grain;;
      506|orientation) run_micro orientation;;
      507|thermo)      run_thermo;;
      508|phase)       run_phase;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_visualization() {
  while true; do
    box
    echo " |                   VISUALIZATION TOOLS                     |"
    box
    echo " | 601) RDF 4x1 merged          604) Elastic constants plot  |"
    echo " | 602) RDF plot                605) Convex hull             |"
    echo " | 603) Born stability 3x2      606) Phase diagram           |"
    echo " | 607) NEP prediction 2x2     609) NEP UMAP+FPS            |"
    echo " | 608) NEP descriptor PCA                                  |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number: "
    read -r opt
    case "$opt" in
      601|rdf4x1)  run_rdf_4x1;;
      602|rdf)     run_rdf_plot;;
      603|bornplot) run_born_plot;;
      604|elasticplot) run_elastic_plot;;
      605|hull)    run_hull;;
      606|exp)     run_phase_exp;;
      607|neppredict) run_nep_predict;;
      608|nepdescriptor) run_nep_descriptor;;
      609|umap)      run_umap;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_utilities() {
  while true; do
    box
    echo " |                      UTILITY TOOLS                       |"
    box
    echo " | 701) QSTEM prepare           704) Show README             |"
    echo " | 702) QSTEM run (multi-core)  705) Structure generator     |"
    echo " | 703) QSTEM list              706) Environment check       |"
    echo " |                             707) Results summary          |"
    box
    echo " | 000) Return to main menu                                    |"
    box
    printf " Input the function number: "
    read -r opt
    case "$opt" in
      701|qprepare) run_qstem_prepare;;
      702|qrun)     run_qstem_run;;
      703|qlist)    run_qstem_list;;
      704|readme)   show_readme;;
      706|env)      run_env_check;;
      707|summary)  run_summary;;
      705|structure) run_structure L12 2;;
      000|0|q|exit) return;;
      *) echo "Unknown: $opt";;
    esac
    pause
  done
}

menu_dev() {
  box
  echo " |                     DEVELOPING / PLANNING                   |"
  box
  echo " |  1) Show README (workflow + roadmap)                       |"
  echo " |  2) List scripts                                            |"
  echo " |  3) Script usage (all modules print usage on error)        |"
  box
  printf " Input the function number: "
  read -r opt
  case "$opt" in
    1|readme) show_readme;;
    2|list)   ls -1 "${SCRIPT_DIR}";;
    3|usage)  echo "Each script prints its Usage when arguments are wrong.";;
    *) echo "Unknown: $opt";;
  esac
  pause
}

# ---------------- main ----------------
main_menu() {
  while true; do
    banner
    echo " ---------------------- HQmdstemkit ----------------------"
    echo " 1) Format Conversion          2) Sample Structures"
    echo " 3) Workflow                   4) Calculators"
    echo " 5) Analyzer                   6) Visualization"
    echo " 7) Utilities                  8) Developing..."
    echo " 0) Exit"
    echo " ------------>>"
    printf " Input the function number: "
    read -r opt || break
    case "$opt" in
      1) menu_format;;
      2) menu_structures;;
      3) menu_workflow;;
      4) menu_calculators;;
      5) menu_analyzer;;
      6) menu_visualization;;
      7) menu_utilities;;
      8) menu_dev;;
      0|q|exit) break;;
      *) echo "Unknown: $opt";;
    esac
  done
}

# direct CLI mode: HQmdstemkit.sh <category> <function> [args...]
if [ $# -ge 2 ]; then
  cat_opt="$1"; fun_opt="$2"; shift 2
  case "${cat_opt}-${fun_opt}" in
    1-101|1-xyz2cfg) run_xyz2cfg "$@";;
    1-102|1-cfg2xyz) run_cfg2xyz "$@";;
    1-103|1-shift)   run_shift "$@";;
    1-104|1-remove)  run_remove "$@";;
    2-201|2-l12)     run_structure L12 "$@";;
    2-202|2-b2)      run_structure B2 "$@";;
    2-203|2-l10)     run_structure L10 "$@";;
    2-204|2-cu80zn20) run_disordered 80 20 "$@";;
    2-205|2-disordered) run_disordered "$@";;
    3-300|3-gpumdsample) run_gpumd_sample "$@";;
    3-301|3-auto)    run_gpumd_auto "$@";;
    3-307|3-prepare) run_gpumd_prepare "$@";;
    3-302|3-run)     run_gpumd_run "$@";;
    3-303|3-train)   run_nep_train "$@";;
    3-304|3-finetune) run_nep_finetune "$@";;
    3-305|3-stem)    run_qstem_stem "$@";;
    3-306|3-elastic0k) run_elastic_0k "$@";;
    3-308|3-batch)     run_batch "$@";;
    3-309|3-abacus)    run_abacus_pretreat "$@";;
    3-310|3-abacusextract) run_abacus_extract "$@";;
    3-311|3-abacusshift)   run_abacus_shift "$@";;
    3-312|3-cp2k)          run_cp2k_pretreat "$@";;
    3-313|3-cp2kextract)   run_cp2k_extract "$@";;
    4-401|4-auto)    run_elastic_auto "$@";;
    4-407|4-strain)  run_elastic_strain "$@";;
    4-402|4-elastic0k) run_elastic_0k "$@";;
    4-403|4-hull)    run_hull "$@";;
    4-404|4-hull-dft) run_hull_dft "$@";;
    4-405|4-exp)     run_phase_exp "$@";;
    4-406|4-born)    run_elastic_born "$@";;
    4-408|4-cohesive) run_cohesive "$@";;
    4-409|4-shear)    run_shear "$@";;
    4-410|4-stackingfault) run_sf "$@";;
    5-501|5-chain)   run_micro chain "$@";;
    5-502|5-cluster) run_micro cluster "$@";;
    5-503|5-segregate) run_micro segregate "$@";;
    5-504|5-twin)    run_micro twin "$@";;
    5-505|5-grain)   run_micro grain "$@";;
    5-506|5-orientation) run_micro orientation "$@";;
    5-507|5-thermo)  run_thermo "$@";;
    5-508|5-phase)   run_phase "$@";;
    6-601|6-rdf4x1)  run_rdf_4x1 "$@";;
    6-602|6-rdf)     run_rdf_plot "$@";;
    6-603|6-bornplot) run_born_plot "$@";;
    6-604|6-elasticplot) run_elastic_plot "$@";;
    6-605|6-hull)    run_hull "$@";;
    6-606|6-exp)     run_phase_exp "$@";;
    6-607|6-neppredict) run_nep_predict "$@";;
    6-608|6-nepdescriptor) run_nep_descriptor "$@";;
    6-609|6-umap)    run_umap "$@";;
    7-701|7-qprepare) run_qstem_prepare "$@";;
    7-702|7-qrun)    run_qstem_run "$@";;
    7-703|7-qlist)   run_qstem_list "$@";;
    7-706|7-env)     run_env_check "$@";;
    7-707|7-summary) run_summary "$@";;
    7-704|7-readme)  show_readme;;
    *) echo "Unknown category/function: ${cat_opt} ${fun_opt}"; echo "Run HQmdstemkit.sh without arguments for the menu."; exit 1;;
  esac
else
  main_menu
fi