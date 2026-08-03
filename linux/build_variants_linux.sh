#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
project_root=$(cd "$linux_dir/.." && pwd)
workspace_root=$(cd "$project_root/.." && pwd)

accelnet_source=${ACCELNET_SOURCE:-$workspace_root/AccelNet}
aenet_source=${AENET_SOURCE:-$workspace_root/aenet-master}
aenet_lammps_source=${AENET_LAMMPS_SOURCE:-$workspace_root/aenet-lammps}
lammps_source=${LAMMPS_SOURCE:-$workspace_root/lammps-4Feb20}
build_jobs=${BENCH_BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}
fortran_compiler=${BENCH_FORTRAN_COMPILER:-gfortran}
run_tests=${BENCH_RUN_TESTS:-1}

require_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $1" >&2
        exit 2
    }
}

require_directory() {
    test -d "$1" || {
        echo "ERROR: required directory not found: $1" >&2
        exit 2
    }
}

for command_name in cmake make mpicxx "$fortran_compiler"; do
    require_command "$command_name"
done
for source_dir in "$accelnet_source" "$aenet_source" "$aenet_lammps_source" "$lammps_source"; do
    require_directory "$source_dir"
done

echo "==> Building original aenet (Release, static, OpenBLAS)"
cmake -S "$aenet_source" -B "$aenet_source/build-linux" \
    -DBUILD_AENET=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
    -DUSE_MPI=OFF \
    -DUSE_OPENBLAS=ON
cmake --build "$aenet_source/build-linux" --target lib --parallel "$build_jobs"

echo "==> Installing aenet and LAMMPS pair-style sources"
install -d "$lammps_source/lib/aenet/include" "$lammps_source/lib/aenet/lib"
install -m 0644 "$aenet_source/src/aenet.h" "$lammps_source/lib/aenet/include/aenet.h"
install -m 0644 "$aenet_source/build-linux/lib/libaenet.a" \
    "$lammps_source/lib/aenet/lib/libaenet.a"
install -m 0644 "$aenet_source/build-linux/lib/liblbfgsb.a" \
    "$lammps_source/lib/aenet/lib/liblbfgsb.a"

install -d "$lammps_source/src/USER-AENET" "$lammps_source/src/USER-ACCELNET"
install -m 0644 "$aenet_lammps_source/USER-AENET/pair_aenet.cpp" \
    "$lammps_source/src/USER-AENET/pair_aenet.cpp"
install -m 0644 "$aenet_lammps_source/USER-AENET/pair_aenet.h" \
    "$lammps_source/src/USER-AENET/pair_aenet.h"
install -m 0755 "$aenet_lammps_source/USER-AENET/Install.sh-gfortran_openblas_serial" \
    "$lammps_source/src/USER-AENET/Install.sh"
install -m 0644 "$project_root/USER-ACCELNET/pair_accelnet.cpp" \
    "$lammps_source/src/USER-ACCELNET/pair_accelnet.cpp"
install -m 0644 "$project_root/USER-ACCELNET/pair_accelnet.h" \
    "$lammps_source/src/USER-ACCELNET/pair_accelnet.h"
install -m 0755 "$project_root/USER-ACCELNET/Install.sh" \
    "$lammps_source/src/USER-ACCELNET/Install.sh"

cd "$lammps_source/src"
install -m 0644 Makefile.package.empty Makefile.package
install -m 0644 Makefile.package.settings.empty Makefile.package.settings
make yes-user-aenet
make yes-user-accelnet
make yes-kspace
make yes-molecule
make yes-rigid

build_runtime_modes() {
    local build_dir="$accelnet_source/build-linux-modes"
    local output_binary="$lammps_source/src/lmp_mpi_accelnet_modes"

    echo "==> Building AccelNet with runtime-selectable auto/direct/moment modes"
    cmake -S "$accelnet_source" -B "$build_dir" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_TESTING=ON \
        -DACCELNET_BUILD_REFERENCE_TESTS=OFF
    cmake --build "$build_dir" --target AccelNetPredictor --parallel "$build_jobs"
    if test "$run_tests" = 1; then
        cmake --build "$build_dir" \
            --target test_descriptor test_tio2 test_lj_model test_behler \
                     test_combined_model test_setup_parser test_n2p2_network \
            --parallel "$build_jobs"
        ctest --test-dir "$build_dir" --output-on-failure \
            -R '^(descriptor_unit|descriptor_tio2|descriptor_lj_model|descriptor_behler|descriptor_combined_model|descriptor_setup_parser|n2p2_model_loader)$'
    fi

    install -d "$lammps_source/lib/accelnet/include" "$lammps_source/lib/accelnet/lib"
    install -m 0644 "$accelnet_source/AccelNetPredictor/include/accelnet.h" \
        "$lammps_source/lib/accelnet/include/accelnet.h"
    install -m 0644 "$build_dir/lib/libaccelnet.a" \
        "$lammps_source/lib/accelnet/lib/libaccelnet.a"
    install -m 0644 "$build_dir/lib/libAccelNetDescriptors.a" \
        "$lammps_source/lib/accelnet/lib/libAccelNetDescriptors.a"

    echo "==> Linking runtime-selectable LAMMPS"
    cd "$lammps_source/src"
    # Recreate the object directory so a fastdep.exe transferred from macOS
    # cannot be reused on Linux (it otherwise fails with Exec format error).
    rm -rf "$lammps_source/src/Obj_mpi"
    rm -f "$lammps_source/src/lmp_mpi"
    make -j"$build_jobs" mpi
    test -x "$lammps_source/src/lmp_mpi" || {
        echo "ERROR: LAMMPS binary was not created" >&2
        exit 3
    }
    install -m 0755 "$lammps_source/src/lmp_mpi" "$output_binary"
    local help_output
    help_output=$("$output_binary" -help 2>&1 || true)
    printf '%s\n' "$help_output" | grep -Eq '(^|[[:space:]])accelnet([[:space:]]|$)' || {
        echo "ERROR: pair_style accelnet missing from $output_binary" >&2
        exit 3
    }
    printf '%s\n' "$help_output" | grep -Eq '(^|[[:space:]])aenet([[:space:]]|$)' || {
        echo "ERROR: pair_style aenet missing from $output_binary" >&2
        exit 3
    }
}

build_runtime_modes

echo
echo "Variant build complete:"
echo "  $lammps_source/src/lmp_mpi_accelnet_modes"
echo "Next: $linux_dir/smoke_test_variants.sh"
