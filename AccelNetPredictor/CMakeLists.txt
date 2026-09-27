cmake_minimum_required(VERSION 3.20)
project(AccelNetPredictor VERSION 1.0.1 LANGUAGES Fortran C)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

include(CMakePackageConfigHelpers)

option(ACCELNET_PREDICTOR_BUILD_TESTING "Build predictor tests" ON)
set(ACCELNET_DESCRIPTORS_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetDescriptors"
    CACHE PATH "AccelNetDescriptors source directory for a same-compiler build")
set(ACCELNET_ORIGINAL_BUILD_DIR "" CACHE PATH
    "Optional out-of-tree AccelNet build used only for the original-speed benchmark")
set(AENET_BUILD_DIR "" CACHE PATH
    "Optional ænet build used only for the direct API benchmark")
set(AENET_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../aenet-master"
    CACHE PATH "Path to the aenet source tree used by reference tests")
set(N2P2_SCALING_EXECUTABLE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../n2p2-master/bin/nnp-scaling"
    CACHE FILEPATH "Path to the n2p2 nnp-scaling reference executable")

if(TARGET AccelNetDescriptors::AccelNetDescriptors)
    # Supplied by a parent superbuild.
elseif(EXISTS "${ACCELNET_DESCRIPTORS_SOURCE_DIR}/CMakeLists.txt")
    set(ACCELNET_DESCRIPTORS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(ACCELNET_DESCRIPTORS_BUILD_REFERENCE OFF CACHE BOOL "" FORCE)
    add_subdirectory("${ACCELNET_DESCRIPTORS_SOURCE_DIR}"
                     "${CMAKE_CURRENT_BINARY_DIR}/AccelNetDescriptors-build")
else()
    find_package(AccelNetDescriptors CONFIG REQUIRED)
endif()
if(AENET_BUILD_DIR AND EXISTS "${AENET_BUILD_DIR}/lib/libaenet.a")
    find_package(BLAS REQUIRED)
    add_executable(aenet-predict-benchmark benchmark/aenet_benchmark.f90)
    target_include_directories(aenet-predict-benchmark PRIVATE
        "${AENET_BUILD_DIR}/modules")
    target_link_libraries(aenet-predict-benchmark PRIVATE
        "${AENET_BUILD_DIR}/lib/libaenet.a" BLAS::BLAS)
    target_compile_options(aenet-predict-benchmark PRIVATE
        $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
    add_executable(aenet-corpus-benchmark benchmark/aenet_corpus_benchmark.f90)
    target_include_directories(aenet-corpus-benchmark PRIVATE
        "${AENET_BUILD_DIR}/modules")
    target_link_libraries(aenet-corpus-benchmark PRIVATE
        "${AENET_BUILD_DIR}/lib/libaenet.a" BLAS::BLAS)
    target_compile_options(aenet-corpus-benchmark PRIVATE
        $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
endif()

add_library(AccelNetPredictor
    src/aenet_network.f90
    src/n2p2_network.f90
    src/predict_input.f90
    src/accelnet_predictor.f90
    src/accelnet_predictor_evaluation.f90
    ../legacy/cpu-reference/accelnet_cpu_reference.f90
    src/accelnet_batch.f90
    src/accelnet.f90)
add_library(AccelNetPredictor::AccelNetPredictor ALIAS AccelNetPredictor)
add_library(AccelNet::AccelNet ALIAS AccelNetPredictor)
set_target_properties(AccelNetPredictor PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME accelnet
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(AccelNetPredictor PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetDescriptors/src/shared"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shared"
    "${CMAKE_CURRENT_SOURCE_DIR}/../legacy/cpu-reference")
target_include_directories(AccelNetPredictor PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(AccelNetPredictor PUBLIC
    AccelNetDescriptors::AccelNetDescriptors)
target_compile_options(AccelNetPredictor PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
    $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)

# Compile the same numerical kernels for the ordinary CPU batch API. Separate
# module names let this serial instance coexist with a GPU-enabled instance.
# Strip directives even when the caller supplies global OpenMP compiler flags.
foreach(_body IN ITEMS generic_force_edge_body composite_batch)
    set(_body_path "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/${_body}.inc")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_body_path}")
    file(READ "${_body_path}" _serial_source)
    string(REGEX REPLACE "[ \t]*!\\$omp[^\n]*" "" _serial_source "${_serial_source}")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/shared")
    configure_file(cmake/serial-source.F90.in
        "${CMAKE_CURRENT_BINARY_DIR}/shared/${_body}_noatomics.inc" @ONLY)
endforeach()
target_include_directories(AccelNetPredictor PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/shared")
set(_shared_modules accelnet_target_runtime accelnet_target_descriptors
    accelnet_target_math accelnet_target_kernels accelnet_batch_target)
foreach(_module IN LISTS _shared_modules)
    set(_extension f90)
    if(_module STREQUAL "accelnet_target_runtime")
        set(_extension F90)
    endif()
    set(_source "${CMAKE_CURRENT_SOURCE_DIR}/src/${_module}.${_extension}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_source}")
    file(READ "${_source}" _serial_source)
    string(REPLACE "'generic_force_edge_body.inc'" "'generic_force_edge_body_noatomics.inc'"
        _serial_source "${_serial_source}")
    string(REPLACE "'composite_batch.inc'" "'composite_batch_noatomics.inc'"
        _serial_source "${_serial_source}")
    foreach(_name IN LISTS _shared_modules)
        string(REPLACE "${_name}" "${_name}_serial" _serial_source "${_serial_source}")
    endforeach()
    string(REGEX REPLACE "[ \t]*!\\$omp[^\n]*" "" _serial_source "${_serial_source}")
    set(_generated "${CMAKE_CURRENT_BINARY_DIR}/serial/${_module}.F90")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/serial")
    # configure_file avoids touching unchanged files and rebuilding everything.
    set(_serial_source "#define ACCELNET_TARGET_SERIAL\n${_serial_source}")
    configure_file(cmake/serial-source.F90.in "${_generated}" @ONLY)
    target_sources(AccelNetPredictor PRIVATE "${_generated}")
endforeach()

option(ACCELNET_BUILD_OPENMP_TARGET "Build the optional OpenMP target batch backend" OFF)
set(ACCELNET_OPENMP_TARGET_FLAGS "" CACHE STRING "Compiler/link flags for OpenMP GPU offload")
option(ACCELNET_TARGET_SERIAL "Compile the target kernels without OpenMP flags/directives for CPU comparison" OFF)
if(ACCELNET_BUILD_OPENMP_TARGET)
    if(NOT ACCELNET_OPENMP_TARGET_FLAGS AND NOT ACCELNET_TARGET_SERIAL)
        message(FATAL_ERROR "Set ACCELNET_OPENMP_TARGET_FLAGS for your compiler/device, e.g. -mp=gpu -gpu=cc90")
    endif()
    set(_offload_flags "")
    if(NOT ACCELNET_TARGET_SERIAL)
        separate_arguments(_offload_flags NATIVE_COMMAND "${ACCELNET_OPENMP_TARGET_FLAGS}")
    endif()
    add_library(AccelNetTarget
        src/accelnet_target_runtime.F90
        src/accelnet_target_descriptors.f90
        src/accelnet_target_math.f90
        src/accelnet_target_kernels.f90
        src/accelnet_batch_target.f90
        src/accelnet_target_c.f90)
    add_library(AccelNet::Target ALIAS AccelNetTarget)
    set_target_properties(AccelNetTarget PROPERTIES OUTPUT_NAME accelnet_target
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/target-modules")
    target_include_directories(AccelNetTarget PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetDescriptors/src/shared"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/shared"
        "${CMAKE_CURRENT_BINARY_DIR}/shared")
    target_include_directories(AccelNetTarget PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/target-modules>"
        "$<INSTALL_INTERFACE:include>")
    target_link_libraries(AccelNetTarget PUBLIC AccelNetPredictor)
    if(ACCELNET_TARGET_SERIAL)
        target_compile_definitions(AccelNetTarget PRIVATE ACCELNET_TARGET_SERIAL)
    endif()
    target_compile_options(AccelNetTarget PRIVATE ${_offload_flags})
    target_link_options(AccelNetTarget PUBLIC ${_offload_flags})
endif()

add_executable(accelnet-predict app/accelnet_predict.f90)
target_link_libraries(accelnet-predict PRIVATE AccelNetPredictor)
add_executable(accelnet-predict-benchmark benchmark/predict_benchmark.f90)
target_link_libraries(accelnet-predict-benchmark PRIVATE AccelNetPredictor)
add_executable(accelnet-corpus-benchmark benchmark/predict_corpus_benchmark.f90)
target_link_libraries(accelnet-corpus-benchmark PRIVATE AccelNetPredictor)
add_executable(accelnet-order-scaling-benchmark benchmark/order_scaling_benchmark.f90)
target_link_libraries(accelnet-order-scaling-benchmark PRIVATE AccelNetPredictor)
if(TARGET aenet-corpus-benchmark)
    add_executable(accelnet-atomic-api-benchmark benchmark/atomic_api_benchmark.f90)
    target_include_directories(accelnet-atomic-api-benchmark PRIVATE
        "${AENET_BUILD_DIR}/modules")
    target_link_libraries(accelnet-atomic-api-benchmark PRIVATE
        AccelNetPredictor "${AENET_BUILD_DIR}/lib/libaenet.a" BLAS::BLAS)
    target_compile_options(accelnet-atomic-api-benchmark PRIVATE
        $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
    add_executable(chebyshev-corpus-validate benchmark/chebyshev_corpus_validate.f90)
    target_include_directories(chebyshev-corpus-validate PRIVATE
        "${AENET_BUILD_DIR}/modules")
    target_link_libraries(chebyshev-corpus-validate PRIVATE
        AccelNetPredictor "${AENET_BUILD_DIR}/lib/libaenet.a" BLAS::BLAS)
    target_compile_options(chebyshev-corpus-validate PRIVATE
        $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
endif()
if(ACCELNET_ORIGINAL_BUILD_DIR AND
   EXISTS "${ACCELNET_ORIGINAL_BUILD_DIR}/lib/libAccelNet-CPU.a")
    add_executable(accelnet-original-predict-benchmark
        benchmark/original_accelnet_benchmark.f90)
    target_include_directories(accelnet-original-predict-benchmark PRIVATE
        "${ACCELNET_ORIGINAL_BUILD_DIR}/mod")
    target_link_libraries(accelnet-original-predict-benchmark PRIVATE
        "${ACCELNET_ORIGINAL_BUILD_DIR}/lib/libAccelNet-CPU.a"
        "${ACCELNET_ORIGINAL_BUILD_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sphericart_c_interface${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set_target_properties(accelnet-original-predict-benchmark PROPERTIES
        BUILD_RPATH "${ACCELNET_ORIGINAL_BUILD_DIR}/lib")
    target_compile_options(accelnet-original-predict-benchmark PRIVATE
        $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
endif()

if(ACCELNET_PREDICTOR_BUILD_TESTING)
    enable_testing()
    option(ACCELNET_TEST_OPENMP_TARGET "Build an opt-in GPU runtime smoke test (not GPU inference)" OFF)
    set(ACCELNET_OPENMP_TARGET_TEST_FLAGS "" CACHE STRING "Compiler/link flags for the GPU runtime test")
    if(ACCELNET_TEST_OPENMP_TARGET)
        if(NOT ACCELNET_OPENMP_TARGET_TEST_FLAGS)
            message(FATAL_ERROR "Set ACCELNET_OPENMP_TARGET_TEST_FLAGS for your compiler and GPU")
        endif()
        separate_arguments(_target_test_flags NATIVE_COMMAND "${ACCELNET_OPENMP_TARGET_TEST_FLAGS}")
        add_executable(test_openmp_target test/test_openmp_target.f90)
        target_compile_options(test_openmp_target PRIVATE ${_target_test_flags})
        target_link_options(test_openmp_target PRIVATE ${_target_test_flags})
        add_test(NAME openmp_target_device_smoke COMMAND test_openmp_target)
        set_tests_properties(openmp_target_device_smoke PROPERTIES
            LABELS "gpu" ENVIRONMENT "OMP_TARGET_OFFLOAD=MANDATORY" TIMEOUT 60)
    endif()
    add_library(batch_test_support STATIC test/batch_test_support.f90)
    set_target_properties(batch_test_support PROPERTIES
        Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/batch-test-modules")
    target_include_directories(batch_test_support PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/batch-test-modules")
    target_link_libraries(batch_test_support PUBLIC AccelNetPredictor)
    if(ACCELNET_BUILD_OPENMP_TARGET)
        add_executable(test_target_c test/test_target_c.c)
        target_include_directories(test_target_c PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
        target_link_libraries(test_target_c PRIVATE AccelNetTarget)
        set_target_properties(test_target_c PROPERTIES LINKER_LANGUAGE Fortran)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
            target_link_options(test_target_c PRIVATE -Mnomain)
        endif()
        if(EXISTS "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii")
            add_test(NAME predictor_target_c_api COMMAND test_target_c
                "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
                "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii")
            set_tests_properties(predictor_target_c_api PROPERTIES LABELS "gpu" TIMEOUT 120)
        endif()
        add_executable(test_batch_target test/test_batch_target.f90)
        target_link_libraries(test_batch_target PRIVATE batch_test_support AccelNetTarget)
        add_executable(accelnet-target-benchmark benchmark/target_benchmark.f90)
        target_link_libraries(accelnet-target-benchmark PRIVATE batch_test_support AccelNetTarget)
        add_test(NAME predictor_target_equivalence COMMAND test_batch_target)
        set_tests_properties(predictor_target_equivalence PROPERTIES LABELS "gpu" TIMEOUT 300)
        foreach(_gpu_mode IN ITEMS direct moment)
            add_test(NAME predictor_target_${_gpu_mode}_equivalence COMMAND test_batch_target --${_gpu_mode})
            set_tests_properties(predictor_target_${_gpu_mode}_equivalence PROPERTIES LABELS "gpu" TIMEOUT 300)
        endforeach()
        add_executable(accelnet-n2p2-benchmark benchmark/n2p2_candidate_benchmark.f90)
        target_link_libraries(accelnet-n2p2-benchmark PRIVATE AccelNetTarget)
        add_executable(write_target_fixtures test/write_target_fixtures.f90)
        target_link_libraries(write_target_fixtures PRIVATE batch_test_support)
        set(_target_fixtures "${CMAKE_CURRENT_BINARY_DIR}/target-fixtures")
        file(MAKE_DIRECTORY "${_target_fixtures}")
        add_test(NAME predictor_target_write_fixtures COMMAND write_target_fixtures "${_target_fixtures}")
        set_tests_properties(predictor_target_write_fixtures PROPERTIES FIXTURES_SETUP target-models LABELS "gpu")
        foreach(_family IN ITEMS lj behler)
            add_test(NAME predictor_target_embedded_${_family} COMMAND test_batch_target --embedded
                "${_target_fixtures}/H.${_family}.nn.ascii" "${_target_fixtures}/O.${_family}.nn.ascii")
            add_test(NAME predictor_target_c_${_family} COMMAND test_target_c
                "${_target_fixtures}/H.${_family}.nn.ascii" "${_target_fixtures}/O.${_family}.nn.ascii")
            set_tests_properties(predictor_target_embedded_${_family} predictor_target_c_${_family}
                PROPERTIES FIXTURES_REQUIRED target-models LABELS "gpu" TIMEOUT 120)
        endforeach()
        add_test(NAME predictor_target_composite COMMAND test_batch_target --composite)
        set_tests_properties(predictor_target_composite PROPERTIES LABELS "gpu" TIMEOUT 300)
        add_test(NAME predictor_target_host_composite COMMAND test_batch_target --host --composite)
        set_tests_properties(predictor_target_host_composite PROPERTIES LABELS "target-host"
            ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
        add_test(NAME predictor_target_g5_moments COMMAND test_batch_target --g5-moments)
        set_tests_properties(predictor_target_g5_moments PROPERTIES LABELS "gpu" TIMEOUT 300)
        add_test(NAME predictor_target_host_g5_moments COMMAND test_batch_target --host --g5-moments)
        set_tests_properties(predictor_target_host_g5_moments PROPERTIES
            LABELS "target-host" ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
        add_test(NAME predictor_target_other_descriptors COMMAND test_batch_target --descriptors)
        set_tests_properties(predictor_target_other_descriptors PROPERTIES LABELS "gpu" TIMEOUT 300)
        add_test(NAME predictor_target_host_equivalence COMMAND test_batch_target --host)
        add_test(NAME predictor_target_host_other_descriptors COMMAND test_batch_target --host --descriptors)
        set_tests_properties(predictor_target_host_equivalence predictor_target_host_other_descriptors
            PROPERTIES LABELS "target-host" ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 300)
        foreach(_n2p2_fixture IN ITEMS n2p2 n2p2-per-element n2p2-per-element-depth n2p2-virial-angular)
            add_test(NAME predictor_target_${_n2p2_fixture} COMMAND test_batch_target --n2p2
                "${CMAKE_CURRENT_SOURCE_DIR}/test/data/${_n2p2_fixture}")
            set_tests_properties(predictor_target_${_n2p2_fixture} PROPERTIES LABELS "gpu" TIMEOUT 300)
        endforeach()
        foreach(_bad_target IN ITEMS offset monotonic target center species shape unsupported activation version10 composite10 uninitialized)
            add_test(NAME predictor_target_reject_${_bad_target} COMMAND test_batch_target ${_bad_target})
            set_tests_properties(predictor_target_reject_${_bad_target} PROPERTIES WILL_FAIL TRUE LABELS "gpu")
        endforeach()
        add_test(NAME predictor_target_reject_host_fallback COMMAND test_batch_target fallback)
        set_tests_properties(predictor_target_reject_host_fallback PROPERTIES
            WILL_FAIL TRUE LABELS "gpu" ENVIRONMENT "OMP_TARGET_OFFLOAD=DISABLED")
    endif()
    add_executable(accelnet-public-api-benchmark benchmark/public_api_benchmark.f90)
    target_link_libraries(accelnet-public-api-benchmark PRIVATE batch_test_support)
    add_executable(accelnet-cpu-regression-benchmark benchmark/cpu_regression_benchmark.F90)
    target_compile_definitions(accelnet-cpu-regression-benchmark PRIVATE ACCELNET_HAVE_BATCH)
    target_link_libraries(accelnet-cpu-regression-benchmark PRIVATE batch_test_support)
    set(ACCELNET_CPU_BASELINE_EXECUTABLE "" CACHE FILEPATH "Archived CPU regression benchmark executable")
    set(ACCELNET_CPU_MAX_SLOWDOWN "1.10" CACHE STRING "Maximum median paired CPU time ratio")
    option(ACCELNET_TEST_COMMON_CPU_PERFORMANCE
        "Compare default common CPU batches against the retained structure evaluator" OFF)
    if(ACCELNET_TEST_COMMON_CPU_PERFORMANCE)
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
        add_test(NAME predictor_common_cpu_performance
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/compare_cpu_performance.py"
                --baseline $<TARGET_FILE:accelnet-cpu-regression-benchmark>
                --candidate $<TARGET_FILE:accelnet-cpu-regression-benchmark>
                --baseline-mode reference --candidate-mode batch
                --max-slowdown "${ACCELNET_CPU_MAX_SLOWDOWN}"
                --output "${CMAKE_CURRENT_BINARY_DIR}/common-cpu-performance.json")
        set_tests_properties(predictor_common_cpu_performance PROPERTIES
            LABELS "performance" RUN_SERIAL TRUE TIMEOUT 300)
    endif()
    option(ACCELNET_TEST_G5_MOMENT_PERFORMANCE
        "Protect the G5 moment advantage using an OpenMP-disabled CPU build" OFF)
    if(ACCELNET_TEST_G5_MOMENT_PERFORMANCE)
        if(NOT ACCELNET_TARGET_SERIAL OR NOT TARGET accelnet-target-benchmark)
            message(FATAL_ERROR "G5 performance tests need the serial target benchmark and OpenMP disabled")
        endif()
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
        add_test(NAME predictor_g5_moment_performance
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/check_g5_moment_performance.py"
                --executable $<TARGET_FILE:accelnet-target-benchmark>
                --max-slowdown "${ACCELNET_CPU_MAX_SLOWDOWN}"
                --output "${CMAKE_CURRENT_BINARY_DIR}/g5-moment-performance.json")
        set_tests_properties(predictor_g5_moment_performance PROPERTIES
            LABELS "performance" RUN_SERIAL TRUE TIMEOUT 300)
    endif()
    if(ACCELNET_CPU_BASELINE_EXECUTABLE)
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
        add_test(NAME predictor_cpu_performance_regression
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/compare_cpu_performance.py"
                --baseline "${ACCELNET_CPU_BASELINE_EXECUTABLE}"
                --candidate $<TARGET_FILE:accelnet-cpu-regression-benchmark>
                --max-slowdown "${ACCELNET_CPU_MAX_SLOWDOWN}"
                --output "${CMAKE_CURRENT_BINARY_DIR}/cpu-performance.json")
        set_tests_properties(predictor_cpu_performance_regression PROPERTIES
            LABELS "performance" RUN_SERIAL TRUE TIMEOUT 300)
    endif()
    add_executable(test_batch test/test_batch.f90)
    target_link_libraries(test_batch PRIVATE batch_test_support)
    add_test(NAME predictor_batch_equivalence COMMAND test_batch "${CMAKE_CURRENT_SOURCE_DIR}/test/data")
    foreach(_bad_batch IN ITEMS offset monotonic target center species shape)
        add_test(NAME predictor_batch_reject_${_bad_batch} COMMAND test_batch
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data" ${_bad_batch})
        set_tests_properties(predictor_batch_reject_${_bad_batch} PROPERTIES WILL_FAIL TRUE)
    endforeach()
    add_executable(test_virial test/test_virial.f90)
    target_link_libraries(test_virial PRIVATE AccelNetPredictor)
    add_test(NAME predictor_atomic_and_structure_virial COMMAND test_virial
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2")
    add_executable(test_virial_c test/test_virial_c.c)
    target_include_directories(test_virial_c PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(test_virial_c PRIVATE AccelNetPredictor)
    add_test(NAME predictor_c_virial COMMAND test_virial_c
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2")
    add_executable(test_virial_convergence test/test_virial_convergence.f90)
    target_link_libraries(test_virial_convergence PRIVATE AccelNetPredictor)
    add_test(NAME predictor_virial_convergence_n2p2 COMMAND test_virial_convergence
        --n2p2 "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2"
        "${CMAKE_CURRENT_BINARY_DIR}/virial-convergence-n2p2.csv")
    add_executable(test_virial_invariants test/test_virial_invariants.f90)
    target_link_libraries(test_virial_invariants PRIVATE AccelNetPredictor)
    foreach(_virial_model IN ITEMS n2p2-per-element n2p2-virial-angular)
        add_test(NAME predictor_virial_convergence_${_virial_model} COMMAND test_virial_convergence
            --n2p2 "${CMAKE_CURRENT_SOURCE_DIR}/test/data/${_virial_model}"
            "${CMAKE_CURRENT_BINARY_DIR}/virial-convergence-${_virial_model}.csv")
        set(_virial_path_check "")
        if(_virial_model STREQUAL "n2p2-virial-angular")
            set(_virial_path_check "mixed-paths")
        endif()
        add_test(NAME predictor_virial_invariants_${_virial_model} COMMAND test_virial_invariants
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/${_virial_model}" ${_virial_path_check})
    endforeach()
    add_executable(test_n2p2_network test/test_n2p2_network.f90)
    target_link_libraries(test_n2p2_network PRIVATE AccelNetPredictor)
    add_test(NAME n2p2_model_loader COMMAND test_n2p2_network
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2")
    add_executable(test_n2p2_topology test/test_n2p2_topology.f90)
    target_link_libraries(test_n2p2_topology PRIVATE AccelNetPredictor)
    add_test(NAME n2p2_per_element_topology_and_normalize_nodes
        COMMAND test_n2p2_topology
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2-per-element"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2-per-element-depth")
    set_tests_properties(n2p2_per_element_topology_and_normalize_nodes
        PROPERTIES LABELS "integration;n2p2-reference")
    add_executable(test_n2p2_atomic_api test/test_n2p2_atomic_api.f90)
    target_link_libraries(test_n2p2_atomic_api PRIVATE AccelNetPredictor)
    add_test(NAME n2p2_atomic_api_and_input_data COMMAND test_n2p2_atomic_api
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2"
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2/input.data")
    add_test(NAME n2p2_input_data_cli COMMAND accelnet-predict --n2p2-data
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2"
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2/input.data")
    set_tests_properties(n2p2_input_data_cli PROPERTIES
        PASS_REGULAR_EXPRESSION "STRUCTURE.*:1;STRUCTURE.*:2")
    if(EXISTS "${N2P2_SCALING_EXECUTABLE}")
        add_executable(dump_n2p2_descriptors test/dump_n2p2_descriptors.f90)
        target_link_libraries(dump_n2p2_descriptors PRIVATE AccelNetPredictor)
        add_executable(compare_n2p2_descriptors test/compare_n2p2_descriptors.f90)
        foreach(N2P2_CUTOFF_TYPE RANGE 0 8)
            add_test(NAME n2p2_descriptor_cutoff_${N2P2_CUTOFF_TYPE}_match
                COMMAND "${CMAKE_COMMAND}"
                    -DN2P2_SCALING=${N2P2_SCALING_EXECUTABLE}
                    -DACCELNET_DUMP=$<TARGET_FILE:dump_n2p2_descriptors>
                    -DCOMPARATOR=$<TARGET_FILE:compare_n2p2_descriptors>
                    -DCASE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/test/n2p2-descriptors
                    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/n2p2-descriptor-${N2P2_CUTOFF_TYPE}
                    -DCUTOFF_TYPE=${N2P2_CUTOFF_TYPE}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_n2p2_scaling.cmake)
            set_tests_properties(n2p2_descriptor_cutoff_${N2P2_CUTOFF_TYPE}_match
                PROPERTIES LABELS "integration;n2p2-reference" TIMEOUT 60)
        endforeach()
    else()
        message(STATUS
            "n2p2 descriptor reference tests disabled: nnp-scaling not found")
    endif()
    add_executable(test_network_activations test/test_network_activations.f90)
    target_link_libraries(test_network_activations PRIVATE AccelNetPredictor)
    add_test(NAME network_activation_compatibility COMMAND test_network_activations
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/legacy-n2p2-softplus.nn.ascii")
    if(EXISTS "${AENET_SOURCE_DIR}/src/ext/feedforward.f90")
        add_library(aenet_feedforward_reference STATIC
            "${AENET_SOURCE_DIR}/src/ext/feedforward.f90")
        set_target_properties(aenet_feedforward_reference PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/aenet-reference-modules")
        target_include_directories(aenet_feedforward_reference PUBLIC
            "${CMAKE_CURRENT_BINARY_DIR}/aenet-reference-modules")
        add_executable(test_aenet_activation_reference
            test/test_aenet_activation_reference.f90)
        target_link_libraries(test_aenet_activation_reference PRIVATE
            AccelNetPredictor aenet_feedforward_reference)
        add_test(NAME network_aenet_activation_reference_match
            COMMAND test_aenet_activation_reference)
        set_tests_properties(network_aenet_activation_reference_match
            PROPERTIES LABELS "aenet-reference")
    endif()
    if(TARGET n2p2_function_reference)
        enable_language(CXX)
        add_executable(test_n2p2_activation_reference
            test/test_n2p2_activation_reference.f90)
        target_link_libraries(test_n2p2_activation_reference PRIVATE
            AccelNetPredictor n2p2_function_reference)
        set_target_properties(test_n2p2_activation_reference PROPERTIES
            LINKER_LANGUAGE CXX)
        add_test(NAME network_n2p2_activation_reference_match
            COMMAND test_n2p2_activation_reference)
        set_tests_properties(network_n2p2_activation_reference_match
            PROPERTIES LABELS "n2p2-reference")
    endif()
    set(ACCELNET_PREDICTOR_GOLDEN_DIR
        "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetTrainer.jl/test/data/fortran_predict"
        CACHE PATH "Existing immutable NN/XSF golden data")
    if(EXISTS "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii")
        if(ACCELNET_BUILD_OPENMP_TARGET)
            add_test(NAME predictor_target_aenet_golden COMMAND test_batch_target
                --aenet "${ACCELNET_PREDICTOR_GOLDEN_DIR}")
            set_tests_properties(predictor_target_aenet_golden PROPERTIES LABELS "gpu" TIMEOUT 300)
            foreach(_gpu_mode RANGE 1 2)
                add_test(NAME predictor_target_aenet_mode${_gpu_mode} COMMAND test_batch_target
                    --aenet "${ACCELNET_PREDICTOR_GOLDEN_DIR}" ${_gpu_mode})
                set_tests_properties(predictor_target_aenet_mode${_gpu_mode} PROPERTIES LABELS "gpu" TIMEOUT 300)
            endforeach()
        endif()
        add_test(NAME predictor_batch_aenet_golden COMMAND test_batch
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data" --aenet "${ACCELNET_PREDICTOR_GOLDEN_DIR}")
        foreach(_virial_structure IN ITEMS structure0001 structure2935)
            add_test(NAME predictor_virial_convergence_${_virial_structure} COMMAND test_virial_convergence
                --aenet "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
                "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
                "${ACCELNET_PREDICTOR_GOLDEN_DIR}/${_virial_structure}.xsf"
                "${CMAKE_CURRENT_BINARY_DIR}/virial-convergence-${_virial_structure}.csv")
        endforeach()
        add_executable(test_network test/test_network.f90)
        target_link_libraries(test_network PRIVATE AccelNetPredictor)
        add_test(NAME network_ascii_roundtrip COMMAND test_network
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii")
        add_executable(test_binary_network test/test_binary_network.f90)
        target_link_libraries(test_binary_network PRIVATE AccelNetPredictor)
        add_test(NAME network_binary_compatibility COMMAND test_binary_network
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
            "${CMAKE_CURRENT_BINARY_DIR}/Ti.test.nn"
            "${CMAKE_CURRENT_BINARY_DIR}/O.test.nn"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf")
        set_tests_properties(network_binary_compatibility PROPERTIES
            FIXTURES_SETUP predictor_binary_networks)
        add_test(NAME predictor_tio2_energy COMMAND accelnet-predict 2
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.fingerprint.stp"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.fingerprint.stp"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf")
        add_executable(test_predictor test/test_predictor.f90)
        target_link_libraries(test_predictor PRIVATE AccelNetPredictor)
        add_test(NAME predictor_fortran_golden COMMAND test_predictor
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.fingerprint.stp"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.fingerprint.stp"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure2935.xsf"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure2935.forces")
        add_test(NAME accelnet_order_scaling_benchmark_smoke
            COMMAND accelnet-order-scaling-benchmark
                "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf" 5)
        set_tests_properties(accelnet_order_scaling_benchmark_smoke PROPERTIES
            PASS_REGULAR_EXPRESSION "ORDER MOMENTS" LABELS performance)
        add_executable(test_accelnet_api test/test_accelnet_api.f90)
        target_link_libraries(test_accelnet_api PRIVATE AccelNetPredictor)
        add_test(NAME accelnet_fortran_api COMMAND test_accelnet_api
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure2935.xsf")
        add_executable(test_accelnet_c_api test/test_accelnet_c_api.c)
        target_include_directories(test_accelnet_c_api PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include")
        target_link_libraries(test_accelnet_c_api PRIVATE AccelNetPredictor)
        add_test(NAME accelnet_c_api COMMAND test_accelnet_c_api
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
            "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2")
        if(TARGET accelnet-atomic-api-benchmark)
            add_test(NAME accelnet_atomic_api_benchmark_smoke
                COMMAND accelnet-atomic-api-benchmark
                    "${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii"
                    "${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii"
                    "${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf" 5)
            set_tests_properties(accelnet_atomic_api_benchmark_smoke PROPERTIES
                PASS_REGULAR_EXPRESSION "ENERGY_FORCE_SPEEDUP_VS_AENET" LABELS performance)
        endif()
        configure_file(test/predict.in.in "${CMAKE_CURRENT_BINARY_DIR}/predict.in" @ONLY)
        add_test(NAME predictor_original_input COMMAND accelnet-predict
            "${CMAKE_CURRENT_BINARY_DIR}/predict.in")
        set_tests_properties(predictor_original_input PROPERTIES
            PASS_REGULAR_EXPRESSION "TOTAL_ENERGY_EV")
        configure_file(test/predict_default.in.in
            "${CMAKE_CURRENT_BINARY_DIR}/predict-default.in" @ONLY)
        add_test(NAME predictor_default_aenet204_version COMMAND accelnet-predict
            "${CMAKE_CURRENT_BINARY_DIR}/predict-default.in")
        set_tests_properties(predictor_default_aenet204_version PROPERTIES
            PASS_REGULAR_EXPRESSION "TOTAL_ENERGY_EV")
        set(ACCELNET_PREDICTOR_CHEBYSHEV_VERSION 1)
        set(ACCELNET_PREDICTOR_VERSION_SEPARATOR " ")
        configure_file(test/predict_version.in.in
            "${CMAKE_CURRENT_BINARY_DIR}/predict-version1.in" @ONLY)
        add_test(NAME predictor_chebyshev_version1 COMMAND accelnet-predict
            "${CMAKE_CURRENT_BINARY_DIR}/predict-version1.in")
        set_tests_properties(predictor_chebyshev_version1 PROPERTIES
            PASS_REGULAR_EXPRESSION "TOTAL_ENERGY_EV")
        set(ACCELNET_PREDICTOR_CHEBYSHEV_VERSION 10)
        set(ACCELNET_PREDICTOR_VERSION_SEPARATOR "=")
        configure_file(test/predict_version.in.in
            "${CMAKE_CURRENT_BINARY_DIR}/predict-version10.in" @ONLY)
        add_test(NAME predictor_chebyshev_version10 COMMAND accelnet-predict
            "${CMAKE_CURRENT_BINARY_DIR}/predict-version10.in")
        set_tests_properties(predictor_chebyshev_version10 PROPERTIES
            PASS_REGULAR_EXPRESSION "TOTAL_ENERGY_EV")
        if(TARGET accelnet-original-predict-benchmark)
            add_test(NAME predictor_original_performance_gate
                COMMAND "${CMAKE_COMMAND}"
                    -DNEW_EXE=$<TARGET_FILE:accelnet-predict-benchmark>
                    -DORIGINAL_EXE=$<TARGET_FILE:accelnet-original-predict-benchmark>
                    -DTI_SETUP=${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.fingerprint.stp
                    -DO_SETUP=${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.fingerprint.stp
                    -DTI_ASCII=${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii
                    -DO_ASCII=${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii
                    -DTI_BINARY=${CMAKE_CURRENT_BINARY_DIR}/Ti.test.nn
                    -DO_BINARY=${CMAKE_CURRENT_BINARY_DIR}/O.test.nn
                    -DSTRUCTURE=${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf
                    -DREPEATS=200
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_original_performance.cmake)
            set_tests_properties(predictor_original_performance_gate PROPERTIES
                FIXTURES_REQUIRED predictor_binary_networks LABELS performance TIMEOUT 60)
        endif()
        if(TARGET aenet-predict-benchmark)
            add_test(NAME predictor_aenet_performance_gate
                COMMAND "${CMAKE_COMMAND}"
                    -DNEW_EXE=$<TARGET_FILE:accelnet-predict-benchmark>
                    -DAENET_EXE=$<TARGET_FILE:aenet-predict-benchmark>
                    -DTI_SETUP=${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.fingerprint.stp
                    -DO_SETUP=${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.fingerprint.stp
                    -DTI_NETWORK=${ACCELNET_PREDICTOR_GOLDEN_DIR}/Ti.nn.ascii
                    -DO_NETWORK=${ACCELNET_PREDICTOR_GOLDEN_DIR}/O.nn.ascii
                    -DSTRUCTURE=${ACCELNET_PREDICTOR_GOLDEN_DIR}/structure0001.xsf
                    -DREPEATS=200
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_aenet_performance.cmake)
            set_tests_properties(predictor_aenet_performance_gate PROPERTIES
                LABELS performance TIMEOUT 60)
        endif()
    endif()
endif()

install(TARGETS AccelNetPredictor accelnet-predict accelnet-predict-benchmark
    accelnet-order-scaling-benchmark
    EXPORT AccelNetPredictorTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/" DESTINATION include OPTIONAL)
if(ACCELNET_BUILD_OPENMP_TARGET)
    install(TARGETS AccelNetTarget EXPORT AccelNetPredictorTargets ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
    install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/target-modules/" DESTINATION include OPTIONAL)
endif()
install(FILES include/accelnet.h include/accelnet_target.h DESTINATION include)
install(EXPORT AccelNetPredictorTargets NAMESPACE AccelNetPredictor::
    DESTINATION lib/cmake/AccelNetPredictor)
configure_package_config_file(cmake/AccelNetPredictorConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetPredictorConfig.cmake"
    INSTALL_DESTINATION lib/cmake/AccelNetPredictor)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetPredictorConfigVersion.cmake"
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/AccelNetPredictorConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetPredictorConfigVersion.cmake"
    DESTINATION lib/cmake/AccelNetPredictor)

# The serial comparison library has no device runtime. Its explicit host tests
# remain enabled; tests requiring an actual GPU cannot run in this build.
if(ACCELNET_TARGET_SERIAL AND ACCELNET_PREDICTOR_BUILD_TESTING)
    get_property(_target_all_tests DIRECTORY PROPERTY TESTS)
    foreach(_test IN LISTS _target_all_tests)
        get_test_property(${_test} LABELS _labels)
        if("gpu" IN_LIST _labels)
            set_tests_properties(${_test} PROPERTIES DISABLED TRUE)
        endif()
    endforeach()
endif()
