#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$linux_dir/.." && pwd)
workspace_root=$(cd "$project_root/.." && pwd)
lammps_binary=${LAMMPS_BINARY:-$workspace_root/lammps-4Feb20/src/lmp_mpi}
result_root=${BENCH_RESULT_DIR:-$linux_dir/results/run-$(date +%Y%m%d-%H%M%S)}

read -r -a systems <<< "${BENCH_SYSTEMS:-water tio2}"
read -r -a sizes <<< "${BENCH_SIZES:-small}"
read -r -a modes <<< "${BENCH_MODES:-fixed dynamic}"
read -r -a pair_styles <<< "${BENCH_PAIR_STYLES:-aenet accelnet}"
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
                "$project_root/examples/water/H.ann" \
                "$project_root/examples/water/O.ann" \
                "$project_root/examples/tio2/models/Ti.15t-15t.nn" \
                "$project_root/examples/tio2/models/O.15t-15t.nn" \
                "$project_root/examples/tio2/data/tio2-relaxed-6.data"
        fi
        for source_dir in "$project_root" "$workspace_root/AccelNet" "$workspace_root/aenet-master" "$workspace_root/lammps-4Feb20"; do
            if git -C "$source_dir" rev-parse HEAD >/dev/null 2>&1; then
                printf 'GIT %s %s\n' "$source_dir" "$(git -C "$source_dir" rev-parse HEAD)"
            fi
        done
        printf 'BENCH_SYSTEMS=%s\n' "${systems[*]}"
        printf 'BENCH_SIZES=%s\n' "${sizes[*]}"
        printf 'BENCH_MODES=%s\n' "${modes[*]}"
        printf 'BENCH_PAIR_STYLES=%s\n' "${pair_styles[*]}"
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
        model_dir="$project_root/examples/tio2/models"
        model_1=Ti.15t-15t.nn
        model_2=O.15t-15t.nn
    else
        echo "ERROR: unsupported system: $system_name" >&2
        exit 2
    fi

    for size_name in "${sizes[@]}"; do
        read -r rep_x rep_y rep_z atom_count <<< "$(replication "$system_name" "$size_name")"
        for mode_name in "${modes[@]}"; do
            input_file="$linux_dir/inputs/in.$system_name.$mode_name"
            for rank_count in "${mpi_ranks[@]}"; do
                for trial in $(seq 1 "$trials"); do
                    for pair_name in "${pair_styles[@]}"; do
                        log_file="$result_root/logs/$system_name-$size_name-$atom_count-$mode_name-$pair_name-r$rank_count-t$trial.log"
                        echo "==> $system_name $size_name $mode_name $pair_name ranks=$rank_count trial=$trial"
                        (
                            cd "$model_dir"
                            run_lammps "$rank_count" "$lammps_binary" -screen none \
                                -var pair_name "$pair_name" \
                                -var data_file "$data_file" \
                                -var model_1 "$model_1" \
                                -var model_2 "$model_2" \
                                -var rep_x "$rep_x" -var rep_y "$rep_y" -var rep_z "$rep_z" \
                                -var warmup_steps "$warmup_steps" \
                                -var measure_steps "$measure_steps" \
                                -in "$input_file" -log "$log_file"
                        )
                    done
                done
            done
        done
    done
done

python3 "$linux_dir/summarize_benchmarks.py" "$result_root/logs" > "$result_root/summary.tsv"
echo "Benchmark completed: $result_root"
echo "Summary: $result_root/summary.tsv"
