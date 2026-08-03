#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$linux_dir/.." && pwd)
workspace_root=$(cd "$project_root/.." && pwd)
lammps_binary=${LAMMPS_BINARY:-$workspace_root/lammps-4Feb20/src/lmp_mpi}
smoke_root=${BENCH_SMOKE_DIR:-$linux_dir/results/smoke-$(date +%Y%m%d-%H%M%S)}

test -x "$lammps_binary" || {
    echo "ERROR: LAMMPS executable not found: $lammps_binary" >&2
    echo "Run $linux_dir/build_linux.sh first." >&2
    exit 2
}
mkdir -p "$smoke_root"

run_system() {
    local system_name=$1
    local data_file model_dir model_1 model_2 input_file
    if test "$system_name" = water; then
        data_file="$project_root/examples/water/01_Start/firstframe.start"
        model_dir="$project_root/examples/water"
        model_1=H.ann
        model_2=O.ann
    else
        data_file="$project_root/examples/tio2/data/tio2-relaxed-6.data"
        model_dir="$project_root/examples/tio2/models"
        model_1=Ti.15t-15t.nn
        model_2=O.15t-15t.nn
    fi
    input_file="$linux_dir/inputs/in.$system_name.run0"

    for pair_name in aenet accelnet; do
        (
            cd "$model_dir"
            "$lammps_binary" -screen none \
                -var pair_name "$pair_name" \
                -var data_file "$data_file" \
                -var model_1 "$model_1" \
                -var model_2 "$model_2" \
                -var dump_file "$smoke_root/$system_name-$pair_name.dump" \
                -in "$input_file" \
                -log "$smoke_root/$system_name-$pair_name.log"
        )
    done

    echo "==> $system_name"
    python3 "$linux_dir/compare_lammps.py" \
        "$smoke_root/$system_name-aenet.log" "$smoke_root/$system_name-aenet.dump" \
        "$smoke_root/$system_name-accelnet.log" "$smoke_root/$system_name-accelnet.dump"
}

run_system water
run_system tio2
echo "Smoke tests completed: $smoke_root"
