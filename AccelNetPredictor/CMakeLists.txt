cmake_minimum_required(VERSION 3.20)
project(AccelNetPredictor VERSION 0.1.0 LANGUAGES Fortran C)
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
    src/accelnet.f90)
add_library(AccelNetPredictor::AccelNetPredictor ALIAS AccelNetPredictor)
add_library(AccelNet::AccelNet ALIAS AccelNetPredictor)
set_target_properties(AccelNetPredictor PROPERTIES
    OUTPUT_NAME accelnet
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(AccelNetPredictor PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
    "$<INSTALL_INTERFACE:include>")
target_link_libraries(AccelNetPredictor PUBLIC
    AccelNetDescriptors::AccelNetDescriptors)
target_compile_options(AccelNetPredictor PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
    $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)

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
    add_executable(test_n2p2_network test/test_n2p2_network.f90)
    target_link_libraries(test_n2p2_network PRIVATE AccelNetPredictor)
    add_test(NAME n2p2_model_loader COMMAND test_n2p2_network
        "${CMAKE_CURRENT_SOURCE_DIR}/test/data/n2p2")
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
install(FILES include/accelnet.h DESTINATION include)
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
