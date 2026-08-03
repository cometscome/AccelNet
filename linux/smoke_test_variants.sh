#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$linux_dir/.." && pwd)
workspace_root=$(cd "$project_root/.." && pwd)
aenet_source=${AENET_SOURCE:-$workspace_root/aenet-master}
lammps_source=${LAMMPS_SOURCE:-$workspace_root/lammps-4Feb20}
lammps_binary=${LAMMPS_ACCELNET_BINARY:-$lammps_source/src/lmp_mpi_accelnet_modes}
smoke_root=${BENCH_SMOKE_DIR:-$linux_dir/results/variant-smoke-$(date +%Y%m%d-%H%M%S)}

test -x "$lammps_binary" || {
    echo "ERROR: LAMMPS executable not found: $lammps_binary" >&2
    echo "Run $linux_dir/build_variants_linux.sh first." >&2
    exit 2
}
mkdir -p "$smoke_root"

run_system() {
    local system_name=$1
    local data_file model_dir model_1 model_2 input_file implementation pair_name backend_name
    if test "$system_name" = water; then
        data_file="$project_root/examples/water/01_Start/firstframe.start"
        model_dir="$project_root/examples/water"
        model_1=H.ann
        model_2=O.ann
    else
        data_file="$project_root/examples/tio2/data/tio2-relaxed-6.data"
        model_dir="$aenet_source/examples/example-chebyshev/03-predict/set001"
        model_1=Ti.15t-15t.nn
        model_2=O.15t-15t.nn
    fi
    test -f "$model_dir/$model_1" -a -f "$model_dir/$model_2" || {
        echo "ERROR: model files not found in $model_dir" >&2
        exit 2
    }
    input_file="$linux_dir/inputs/in.$system_name.run0"

    for implementation in aenet accelnet-direct accelnet-moment; do
        case "$implementation" in
            aenet) pair_name=aenet; backend_name=auto ;;
            accelnet-direct) pair_name=accelnet; backend_name=direct ;;
            accelnet-moment) pair_name=accelnet; backend_name=moment ;;
        esac
        (
            cd "$model_dir"
            "$lammps_binary" -screen none \
                -var pair_name "$pair_name" \
                -var backend_name "$backend_name" \
                -var data_file "$data_file" \
                -var model_1 "$model_1" \
                -var model_2 "$model_2" \
                -var dump_file "$smoke_root/$system_name-$implementation.dump" \
                -in "$input_file" \
                -log "$smoke_root/$system_name-$implementation.log"
        )
    done

    grep -q '^AccelNet Chebyshev evaluation mode: direct$' \
        "$smoke_root/$system_name-accelnet-direct.log" || {
        echo "ERROR: direct mode was not confirmed in the LAMMPS log" >&2
        exit 3
    }
    grep -q '^AccelNet Chebyshev evaluation mode: moment$' \
        "$smoke_root/$system_name-accelnet-moment.log" || {
        echo "ERROR: moment mode was not confirmed in the LAMMPS log" >&2
        exit 3
    }

    echo "==> $system_name: aenet vs AccelNet direct"
    python3 "$linux_dir/compare_lammps.py" \
        "$smoke_root/$system_name-aenet.log" "$smoke_root/$system_name-aenet.dump" \
        "$smoke_root/$system_name-accelnet-direct.log" "$smoke_root/$system_name-accelnet-direct.dump"
    echo "==> $system_name: aenet vs AccelNet moment"
    python3 "$linux_dir/compare_lammps.py" \
        "$smoke_root/$system_name-aenet.log" "$smoke_root/$system_name-aenet.dump" \
        "$smoke_root/$system_name-accelnet-moment.log" "$smoke_root/$system_name-accelnet-moment.dump"
}

run_system water
run_system tio2
echo "Variant smoke tests completed: $smoke_root"
