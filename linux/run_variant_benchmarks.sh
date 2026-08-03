#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$linux_dir/.." && pwd)
workspace_root=$(cd "$project_root/.." && pwd)
aenet_source=${AENET_SOURCE:-$workspace_root/aenet-master}
lammps_source=${LAMMPS_SOURCE:-$workspace_root/lammps-4Feb20}
lammps_binary=${LAMMPS_ACCELNET_BINARY:-$lammps_source/src/lmp_mpi_accelnet_modes}
result_root=${BENCH_RESULT_DIR:-$linux_dir/results/variant-run-$(date +%Y%m%d-%H%M%S)}

read -r -a systems <<< "${BENCH_SYSTEMS:-water tio2}"
read -r -a sizes <<< "${BENCH_SIZES:-small}"
read -r -a modes <<< "${BENCH_MODES:-fixed dynamic}"
read -r -a implementations <<< "${BENCH_IMPLEMENTATIONS:-aenet accelnet-direct accelnet-moment}"
read -r -a mpi_ranks <<< "${BENCH_MPI_RANKS:-1 2}"
trials=${BENCH_TRIALS:-3}
warmup_steps=${BENCH_WARMUP_STEPS:-100}
measure_steps=${BENCH_MEASURE_STEPS:-500}
launcher=${BENCH_LAUNCHER:-mpirun}
bind_cores=${BENCH_BIND_CORES:-1}

export OMP_NUM_THREADS=${BENCH_OMP_THREADS:-1}
export OPENBLAS_NUM_THREADS=${BENCH_BLAS_THREADS:-1}
export MKL_NUM_THREADS=${BENCH_BLAS_THREADS:-1}
export BLIS_NUM_THREADS=${BENCH_BLAS_THREADS:-1}

test -x "$lammps_binary" || {
    echo "ERROR: LAMMPS executable not found: $lammps_binary" >&2
    exit 2
}
command -v "$launcher" >/dev/null 2>&1 || {
    echo "ERROR: launcher not found: $launcher" >&2
    exit 2
}
mkdir -p "$result_root/logs"
printf 'RUNNING\n' > "$result_root/status.txt"
printf '%s\n' "$$" > "$result_root/runner.pid"

finish_run() {
    local exit_code=$?
    if test "$exit_code" = 0; then
        printf 'COMPLETED\n' > "$result_root/status.txt"
    else
        printf 'FAILED exit_code=%s\n' "$exit_code" > "$result_root/status.txt"
    fi
}
trap finish_run EXIT
interrupt_run() {
    local exit_code=$1
    printf 'INTERRUPTED exit_code=%s\n' "$exit_code" > "$result_root/status.txt"
    trap - EXIT
    exit "$exit_code"
}
trap 'interrupt_run 130' INT
trap 'interrupt_run 143' TERM

replication() {
    case "$1:$2" in
        water:small)  echo "3 2 1 1152" ;;
        water:medium) echo "6 5 2 11520" ;;
        water:large)  echo "9 8 7 96768" ;;
        tio2:small)   echo "8 6 4 1152" ;;
        tio2:medium)  echo "12 14 8 8064" ;;
        tio2:large)   echo "32 32 16 98304" ;;
        *) echo "ERROR: unsupported system/size: $1/$2" >&2; return 2 ;;
    esac
}

select_implementation() {
    case "$1" in
        aenet) selected_pair=aenet; selected_backend=auto ;;
        accelnet-direct) selected_pair=accelnet; selected_backend=direct ;;
        accelnet-moment) selected_pair=accelnet; selected_backend=moment ;;
        *) echo "ERROR: unsupported implementation: $1" >&2; return 2 ;;
    esac
}

write_system_info() {
    {
        date --iso-8601=seconds 2>/dev/null || date
        uname -a
        command -v lscpu >/dev/null 2>&1 && lscpu
        "$lammps_binary" -help 2>&1 | head -40 || true
        mpicxx --version 2>&1 | head -5 || true
        gfortran --version 2>&1 | head -5 || true
        "$launcher" --version 2>&1 | head -10 || true
        cmake --version 2>&1 | head -3 || true
        if command -v sha256sum >/dev/null 2>&1; then
            sha256sum \
                "$lammps_binary" \
                "$workspace_root/AccelNet/AccelNetDescriptors/CMakeLists.txt" \
                "$workspace_root/AccelNet/AccelNetDescriptors/src/accelnet_descriptors.f90" \
                "$project_root/examples/water/H.ann" \
                "$project_root/examples/water/O.ann" \
                "$aenet_source/examples/example-chebyshev/03-predict/set001/Ti.15t-15t.nn" \
                "$aenet_source/examples/example-chebyshev/03-predict/set001/O.15t-15t.nn" \
                "$project_root/examples/tio2/data/tio2-relaxed-6.data"
        fi
        for source_dir in "$project_root" "$workspace_root/AccelNet" "$aenet_source" "$lammps_source"; do
            if git -C "$source_dir" rev-parse HEAD >/dev/null 2>&1; then
                printf 'GIT %s %s\n' "$source_dir" "$(git -C "$source_dir" rev-parse HEAD)"
                if test -n "$(git -C "$source_dir" status --short 2>/dev/null)"; then
                    printf 'GIT_DIRTY %s yes\n' "$source_dir"
                fi
            fi
        done
        printf 'ACCELNET_RUNTIME_MODES=auto direct moment\n'
        printf 'BENCH_SYSTEMS=%s\n' "${systems[*]}"
        printf 'BENCH_SIZES=%s\n' "${sizes[*]}"
        printf 'BENCH_MODES=%s\n' "${modes[*]}"
        printf 'BENCH_IMPLEMENTATIONS=%s\n' "${implementations[*]}"
        printf 'BENCH_MPI_RANKS=%s\n' "${mpi_ranks[*]}"
        printf 'BENCH_TRIALS=%s\n' "$trials"
        printf 'BENCH_WARMUP_STEPS=%s\n' "$warmup_steps"
        printf 'BENCH_MEASURE_STEPS=%s\n' "$measure_steps"
        printf 'BENCH_BIND_CORES=%s\n' "$bind_cores"
        printf 'OMP_NUM_THREADS=%s\n' "$OMP_NUM_THREADS"
        printf 'OPENBLAS_NUM_THREADS=%s\n' "$OPENBLAS_NUM_THREADS"
    } > "$result_root/system-info.txt"
}

run_lammps() {
    local rank_count=$1
    shift
    local -a command_line
    if test "$launcher" = srun; then
        command_line=(srun -n "$rank_count")
        test "$bind_cores" = 1 && command_line+=(--cpu-bind=cores)
    else
        command_line=("$launcher" -np "$rank_count")
        test "$bind_cores" = 1 && command_line+=(--bind-to core --map-by core)
    fi
    command_line+=("$@")
    "${command_line[@]}"
}

write_system_info
for system_name in "${systems[@]}"; do
    if test "$system_name" = water; then
        data_file="$project_root/examples/water/01_Start/firstframe.start"
        model_dir="$project_root/examples/water"
        model_1=H.ann
        model_2=O.ann
    elif test "$system_name" = tio2; then
        data_file="$project_root/examples/tio2/data/tio2-relaxed-6.data"
        model_dir="$aenet_source/examples/example-chebyshev/03-predict/set001"
        model_1=Ti.15t-15t.nn
        model_2=O.15t-15t.nn
    else
        echo "ERROR: unsupported system: $system_name" >&2
        exit 2
    fi
    test -f "$model_dir/$model_1" -a -f "$model_dir/$model_2" || {
        echo "ERROR: model files not found in $model_dir" >&2
        exit 2
    }

    for size_name in "${sizes[@]}"; do
        read -r rep_x rep_y rep_z atom_count <<< "$(replication "$system_name" "$size_name")"
        for mode_name in "${modes[@]}"; do
            input_file="$linux_dir/inputs/in.$system_name.$mode_name"
            for rank_count in "${mpi_ranks[@]}"; do
                for trial in $(seq 1 "$trials"); do
                    for implementation in "${implementations[@]}"; do
                        select_implementation "$implementation"
                        log_file="$result_root/logs/$system_name-$size_name-$atom_count-$mode_name-$implementation-r$rank_count-t$trial.log"
                        echo "==> $system_name $size_name $mode_name $implementation ranks=$rank_count trial=$trial"
                        (
                            cd "$model_dir"
                            run_lammps "$rank_count" "$lammps_binary" -screen none \
                                -var pair_name "$selected_pair" \
                                -var backend_name "$selected_backend" \
                                -var data_file "$data_file" \
                                -var model_1 "$model_1" \
                                -var model_2 "$model_2" \
                                -var rep_x "$rep_x" -var rep_y "$rep_y" -var rep_z "$rep_z" \
                                -var warmup_steps "$warmup_steps" \
                                -var measure_steps "$measure_steps" \
                                -in "$input_file" -log "$log_file"
                        )
                        if test "$selected_pair" = accelnet; then
                            grep -q "^AccelNet Chebyshev evaluation mode: $selected_backend$" \
                                "$log_file" || {
                                echo "ERROR: $selected_backend mode was not confirmed in $log_file" >&2
                                exit 3
                            }
                        fi
                    done
                done
            done
        done
    done
done

python3 "$linux_dir/summarize_variant_benchmarks.py" "$result_root/logs" > "$result_root/summary.tsv"
echo "Benchmark completed: $result_root"
echo "Summary: $result_root/summary.tsv"
