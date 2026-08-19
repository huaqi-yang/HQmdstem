#!/bin/bash
# HQmdstemkit: ABACUS OUT.ABACUS -> total_train.xyz (NEP-ready extxyz)
output_file="total_train.xyz"
> "$output_file"

directories=$(find . -maxdepth 3 -type d -name "OUT.ABACUS" | xargs -I {} dirname {})

count=0
total=$(echo "$directories" | wc -l)

for dir in $directories; do
    count=$((count + 1))
    echo "Processing [$count/$total]: $dir"
    python3 - "$dir" "$output_file" <<'EOF'
import sys, dpdata
from ase.io import write
work_dir, out_file = sys.argv[1], sys.argv[2]
try:
    sys_ = dpdata.LabeledSystem(work_dir, fmt='abacus/scf')
    if len(sys_) > 0:
        traj = sys_.to_ase_structure()
        write(out_file, traj, format='extxyz', append=True)
except Exception as e:
    print(f"skip {work_dir}: {e}")
EOF
done

echo "Done: merged ABACUS frames into $output_file"