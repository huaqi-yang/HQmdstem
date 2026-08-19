PSEUDO_DIR="$HOME/CuZn"
ORB_DIR="$HOME/CuZn"

for i in $(seq 1 383); do
    file="POSCAR_${i}.vasp"
    target_dir="iter00_${i}"
    
    if [ -f "$file" ]; then
        echo "Processing $file -> $target_dir"
        
        abacustest model inputs -f "$file" --ftype poscar --lcao --jtype scf
        
        if [ ! -d "000000" ]; then
            echo "ERROR: abacustest failed for $file"
            continue
        fi
        
        rm -rf "$target_dir"
        mv 000000 "$target_dir"
        
        cat <<EOF > "$target_dir/STRU.tmp"
ATOMIC_SPECIES
Cu 63.546 Cu_ONCV_PBE-1.0.upf
Zn 65.38 Zn_ONCV_PBE_FR-1.0.upf
NUMERICAL_ORBITAL
Cu_gga_10au_150Ry_6s3p3d2f.orb
Zn_gga_10au_150Ry_6s3p3d2f.orb
EOF
        sed -n '/LATTICE_CONSTANT/,$p' "$target_dir/STRU" >> "$target_dir/STRU.tmp"
        mv "$target_dir/STRU.tmp" "$target_dir/STRU"
        
        cat <<EOF > "$target_dir/INPUT"
INPUT_PARAMETERS
pseudo_dir      ${PSEUDO_DIR}
orbital_dir     ${ORB_DIR}
basis_type      lcao
ecutwfc         100
scf_nmax        100
scf_thr         1e-6
device          gpu
ks_solver       cusolver
precision       double
cal_force       1
cal_stress      1
kspacing        0.14
smearing_method gaussian
smearing_sigma  0.02
EOF
    else
        echo "Skip: $file not found"
    fi
done
