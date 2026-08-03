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

echo "==> Building AccelNet (Release, static)"
cmake -S "$accelnet_source" -B "$accelnet_source/build-linux" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
    -DBUILD_SHARED_LIBS=OFF \
    -DBUILD_TESTING=ON \
    -DACCELNET_BUILD_REFERENCE_TESTS=OFF
cmake --build "$accelnet_source/build-linux" \
    --target AccelNetPredictor --parallel "$build_jobs"
if test "$run_tests" = 1; then
    cmake --build "$accelnet_source/build-linux" \
        --target test_descriptor test_tio2 test_lj_model test_behler \
                 test_combined_model test_setup_parser test_n2p2_network \
        --parallel "$build_jobs"
    ctest --test-dir "$accelnet_source/build-linux" --output-on-failure \
        -R '^(descriptor_unit|descriptor_tio2|descriptor_lj_model|descriptor_behler|descriptor_combined_model|descriptor_setup_parser|n2p2_model_loader)$'
fi

echo "==> Building original aenet (Release, static, OpenBLAS)"
cmake -S "$aenet_source" -B "$aenet_source/build-linux" \
    -DBUILD_AENET=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
    -DUSE_MPI=OFF \
    -DUSE_OPENBLAS=ON
cmake --build "$aenet_source/build-linux" --target lib --parallel "$build_jobs"

echo "==> Installing libraries into the LAMMPS source tree"
install -d "$lammps_source/lib/accelnet/include" "$lammps_source/lib/accelnet/lib"
install -m 0644 "$accelnet_source/AccelNetPredictor/include/accelnet.h" \
    "$lammps_source/lib/accelnet/include/accelnet.h"
install -m 0644 "$accelnet_source/build-linux/lib/libaccelnet.a" \
    "$lammps_source/lib/accelnet/lib/libaccelnet.a"
install -m 0644 "$accelnet_source/build-linux/lib/libAccelNetDescriptors.a" \
    "$lammps_source/lib/accelnet/lib/libAccelNetDescriptors.a"

install -d "$lammps_source/lib/aenet/include" "$lammps_source/lib/aenet/lib"
install -m 0644 "$aenet_source/src/aenet.h" "$lammps_source/lib/aenet/include/aenet.h"
install -m 0644 "$aenet_source/build-linux/lib/libaenet.a" \
    "$lammps_source/lib/aenet/lib/libaenet.a"
install -m 0644 "$aenet_source/build-linux/lib/liblbfgsb.a" \
    "$lammps_source/lib/aenet/lib/liblbfgsb.a"

echo "==> Installing USER-AENET and USER-ACCELNET pair styles"
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

echo "==> Building combined LAMMPS executable"
cd "$lammps_source/src"
make yes-user-aenet
make yes-user-accelnet
make yes-kspace
make yes-molecule
make yes-rigid
make -j"$build_jobs" mpi

lammps_binary="$lammps_source/src/lmp_mpi"
test -x "$lammps_binary" || {
    echo "ERROR: LAMMPS binary was not created: $lammps_binary" >&2
    exit 3
}
help_output=$("$lammps_binary" -help 2>&1 || true)
printf '%s\n' "$help_output" | grep -Eq '(^|[[:space:]])accelnet([[:space:]]|$)' || {
    echo "ERROR: pair_style accelnet is missing from the built executable" >&2
    exit 3
}
printf '%s\n' "$help_output" | grep -Eq '(^|[[:space:]])aenet([[:space:]]|$)' || {
    echo "ERROR: pair_style aenet is missing from the built executable" >&2
    exit 3
}

echo
echo "Build complete: $lammps_binary"
echo "Next: $linux_dir/smoke_test.sh"
