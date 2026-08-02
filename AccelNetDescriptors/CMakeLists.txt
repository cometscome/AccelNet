cmake_minimum_required(VERSION 3.20)
project(AccelNetDescriptors VERSION 0.1.0 LANGUAGES Fortran)

include(CMakePackageConfigHelpers)

option(ACCELNET_DESCRIPTORS_BUILD_TESTING "Build tests" ON)
option(ACCELNET_DESCRIPTORS_BUILD_REFERENCE "Build comparisons against an untouched AccelNet tree" ON)
set(ACCELNET_CHEBYSHEV_BACKEND "auto" CACHE STRING
    "Chebyshev angular backend: auto, direct, or moment")
set_property(CACHE ACCELNET_CHEBYSHEV_BACKEND PROPERTY STRINGS auto direct moment)
string(TOLOWER "${ACCELNET_CHEBYSHEV_BACKEND}" ACCELNET_CHEBYSHEV_BACKEND_NORMALIZED)
if(NOT ACCELNET_CHEBYSHEV_BACKEND_NORMALIZED MATCHES "^(auto|direct|moment)$")
    message(FATAL_ERROR
        "ACCELNET_CHEBYSHEV_BACKEND must be auto, direct, or moment")
endif()
set(ACCELNET_ORIGINAL_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNet-a76ec3a83c8b0246328e9cca4817cc0afc84c925"
    CACHE PATH "Path to the untouched AccelNet source tree")

add_library(AccelNetDescriptors
    src/accelnet_legacy_lcl.f90
    src/accelnet_descriptors.f90
    src/accelnet_lj.f90
    src/accelnet_behler.f90
    src/accelnet_descriptor_models.f90
    src/accelnet_setup.f90
)
set_source_files_properties(src/accelnet_descriptors.f90 PROPERTIES
    Fortran_PREPROCESS ON)
add_library(AccelNetDescriptors::AccelNetDescriptors ALIAS AccelNetDescriptors)
set_target_properties(AccelNetDescriptors PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(AccelNetDescriptors PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
    "$<INSTALL_INTERFACE:include>")
target_compile_options(AccelNetDescriptors PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
    $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)
if(ACCELNET_CHEBYSHEV_BACKEND_NORMALIZED STREQUAL "direct")
    target_compile_definitions(AccelNetDescriptors PRIVATE
        ACCELNET_CHEBYSHEV_FORCE_DIRECT=1)
elseif(ACCELNET_CHEBYSHEV_BACKEND_NORMALIZED STREQUAL "moment")
    target_compile_definitions(AccelNetDescriptors PRIVATE
        ACCELNET_CHEBYSHEV_FORCE_MOMENT=1)
endif()

add_executable(accelnet-descriptor app/descriptor_generate.f90)
target_link_libraries(accelnet-descriptor PRIVATE AccelNetDescriptors)
add_executable(accelnet-setup-descriptor app/setup_descriptor_generate.f90)
target_link_libraries(accelnet-setup-descriptor PRIVATE AccelNetDescriptors)
add_executable(benchmark-g4-derivative benchmark/benchmark_g4_derivative.f90)
target_compile_options(benchmark-g4-derivative PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>)
add_executable(benchmark-g5-scaling benchmark/benchmark_g5_scaling.f90)
target_link_libraries(benchmark-g5-scaling PRIVATE AccelNetDescriptors)

if(ACCELNET_DESCRIPTORS_BUILD_TESTING)
    enable_testing()

    add_executable(test_descriptor test/test_descriptor.f90)
    target_link_libraries(test_descriptor PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_unit COMMAND test_descriptor)

    add_executable(test_tio2 test/test_tio2.f90)
    target_link_libraries(test_tio2 PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_tio2
        COMMAND test_tio2
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf")

    add_executable(test_lj_model test/test_lj_model.f90)
    target_link_libraries(test_lj_model PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_lj_model COMMAND test_lj_model)

    add_executable(test_behler test/test_behler.f90)
    target_link_libraries(test_behler PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_behler COMMAND test_behler)

    add_executable(test_combined_model test/test_combined_model.f90)
    target_link_libraries(test_combined_model PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_combined_model COMMAND test_combined_model)

    add_executable(dump_descriptor_reference test/dump_descriptor_reference.f90)
    target_link_libraries(dump_descriptor_reference PRIVATE AccelNetDescriptors)

    add_executable(test_setup_parser test/test_setup_parser.f90)
    target_link_libraries(test_setup_parser PRIVATE AccelNetDescriptors)
    add_test(NAME descriptor_setup_parser
        COMMAND test_setup_parser
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Ti.fingerprint.stp"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/LJ.fingerprint.stp"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Behler.fingerprint.stp"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Multi.fingerprint.stp"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v1.fingerprint.stp"
            "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v10.fingerprint.stp")
    add_test(NAME descriptor_setup_cli_match
        COMMAND "${CMAKE_COMMAND}"
            -DNEW_EXE=$<TARGET_FILE:accelnet-setup-descriptor>
            -DLEGACY_EXE=$<TARGET_FILE:accelnet-descriptor>
            -DTI_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/Ti.fingerprint.stp
            -DO_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/O.fingerprint.stp
            -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf
            -DWORK=${CMAKE_CURRENT_BINARY_DIR}/setup-cli-comparison
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_setup_cli.cmake)

    if(ACCELNET_DESCRIPTORS_BUILD_REFERENCE)
        if(NOT EXISTS "${ACCELNET_ORIGINAL_DIR}/src/chebyshev/fingerprint_basis_chebyshev.f90")
            message(FATAL_ERROR "Untouched AccelNet source not found: ${ACCELNET_ORIGINAL_DIR}")
        endif()

        add_library(accelnet_original_chebyshev STATIC
            "${ACCELNET_ORIGINAL_DIR}/src/mod_fingerprint_basis.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/utils/cutoff.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/chebyshev/chebyshev.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/chebyshev/fingerprint_basis_chebyshev.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/LJ/fingerprint_basis_LJ.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/symmetryfunctions/cutoff_cache.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/symmetryfunctions/sf_indices.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/symmetryfunctions/behler_functions.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/symmetryfunctions/fingerprint_basis_behler2011.f90")
        set_target_properties(accelnet_original_chebyshev PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/original-modules")
        target_include_directories(accelnet_original_chebyshev PUBLIC
            "${CMAKE_CURRENT_BINARY_DIR}/original-modules")
        target_compile_options(accelnet_original_chebyshev PRIVATE
            $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
            $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)

        add_library(accelnet_original_lcl STATIC
            "${ACCELNET_ORIGINAL_DIR}/src/utils/sortlib.f90"
            "${ACCELNET_ORIGINAL_DIR}/src/utils/lclist.f90")
        set_target_properties(accelnet_original_lcl PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/original-lcl-modules")
        target_include_directories(accelnet_original_lcl PUBLIC
            "${CMAKE_CURRENT_BINARY_DIR}/original-lcl-modules")
        target_compile_options(accelnet_original_lcl PRIVATE
            $<$<Fortran_COMPILER_ID:GNU>:-O3;-ffree-line-length-none>
            $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O3>)

        add_executable(original_descriptor_output
            test/original_descriptor_output.f90)
        target_link_libraries(original_descriptor_output PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev accelnet_original_lcl)

        add_executable(compare_descriptor_files
            test/compare_descriptor_files.f90)

        add_executable(test_lj_reference
            test/test_lj_reference.f90)
        target_link_libraries(test_lj_reference PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev)
        add_test(NAME descriptor_lj_original_match
            COMMAND test_lj_reference
                "${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf")

        add_executable(test_chebyshev_versions_reference
            test/test_chebyshev_versions_reference.f90)
        target_link_libraries(test_chebyshev_versions_reference PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev)
        add_test(NAME descriptor_chebyshev_versions_original_match
            COMMAND test_chebyshev_versions_reference)

        add_executable(test_behler_reference
            test/test_behler_reference.f90)
        target_link_libraries(test_behler_reference PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev)
        add_test(NAME descriptor_behler_original_match
            COMMAND test_behler_reference
                "${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf"
                "${CMAKE_CURRENT_SOURCE_DIR}/test/data/Behler.fingerprint.stp")

        add_executable(compare_performance
            benchmark/compare_performance.f90)
        target_link_libraries(compare_performance PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev accelnet_original_lcl)

        add_executable(compare_extension_performance
            benchmark/compare_extension_performance.f90)
        target_link_libraries(compare_extension_performance PRIVATE
            AccelNetDescriptors accelnet_original_chebyshev)

        add_test(NAME descriptor_original_file_match
            COMMAND "${CMAKE_COMMAND}"
                -DNEW_EXE=$<TARGET_FILE:accelnet-descriptor>
                -DORIGINAL_EXE=$<TARGET_FILE:original_descriptor_output>
                -DCOMPARATOR=$<TARGET_FILE:compare_descriptor_files>
                -DINPUTS=${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/file-comparison
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_outputs.cmake)

        add_test(NAME descriptor_chebyshev_v1_file_match
            COMMAND "${CMAKE_COMMAND}"
                -DNEW_EXE=$<TARGET_FILE:accelnet-setup-descriptor>
                -DORIGINAL_EXE=$<TARGET_FILE:original_descriptor_output>
                -DCOMPARATOR=$<TARGET_FILE:compare_descriptor_files>
                -DVERSION=1
                -DTI_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v1.fingerprint.stp
                -DO_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v1-O.fingerprint.stp
                -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/chebyshev-version-file-comparison
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_setup_versions.cmake)
        add_test(NAME descriptor_chebyshev_v10_file_match
            COMMAND "${CMAKE_COMMAND}"
                -DNEW_EXE=$<TARGET_FILE:accelnet-setup-descriptor>
                -DORIGINAL_EXE=$<TARGET_FILE:original_descriptor_output>
                -DCOMPARATOR=$<TARGET_FILE:compare_descriptor_files>
                -DVERSION=10
                -DTI_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v10.fingerprint.stp
                -DO_SETUP=${CMAKE_CURRENT_SOURCE_DIR}/test/data/Chebyshev-v10-O.fingerprint.stp
                -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/test/data/structure0001.xsf
                -DWORK=${CMAKE_CURRENT_BINARY_DIR}/chebyshev-version-file-comparison
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_setup_versions.cmake)

        set(TIO2_DATASET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../TiO2-xsf"
            CACHE PATH "Optional TiO2 XSF corpus for integration tests")
        if(EXISTS "${TIO2_DATASET_DIR}/structure7815.xsf")
            add_test(NAME descriptor_original_tio2_corpus_match
                COMMAND "${CMAKE_COMMAND}"
                    -DNEW_EXE=$<TARGET_FILE:accelnet-descriptor>
                    -DORIGINAL_EXE=$<TARGET_FILE:original_descriptor_output>
                    -DCOMPARATOR=$<TARGET_FILE:compare_descriptor_files>
                    "-DINPUTS=${TIO2_DATASET_DIR}/structure0001.xsf,${TIO2_DATASET_DIR}/structure1000.xsf,${TIO2_DATASET_DIR}/structure3000.xsf,${TIO2_DATASET_DIR}/structure5000.xsf,${TIO2_DATASET_DIR}/structure7000.xsf,${TIO2_DATASET_DIR}/structure7815.xsf"
                    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/tio2-corpus-comparison
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/compare_outputs.cmake)
            add_test(NAME descriptor_performance_gate
                COMMAND compare_performance
                    "${TIO2_DATASET_DIR}/structure7513.xsf" 50)
            add_test(NAME descriptor_chebyshev_v1_performance_gate
                COMMAND compare_performance
                    "${TIO2_DATASET_DIR}/structure7513.xsf" 50 1)
            add_test(NAME descriptor_chebyshev_v10_performance_gate
                COMMAND compare_performance
                    "${TIO2_DATASET_DIR}/structure7513.xsf" 50 10)
            add_test(NAME descriptor_extension_performance_gate
                COMMAND compare_extension_performance
                    "${TIO2_DATASET_DIR}/structure7513.xsf" 200)
            set_tests_properties(descriptor_performance_gate PROPERTIES
                LABELS performance
                TIMEOUT 120)
            set_tests_properties(
                descriptor_chebyshev_v1_performance_gate
                descriptor_chebyshev_v10_performance_gate
                PROPERTIES LABELS performance TIMEOUT 120)
            set_tests_properties(descriptor_extension_performance_gate PROPERTIES
                LABELS performance
                TIMEOUT 120)
        endif()
    endif()
endif()

install(TARGETS AccelNetDescriptors accelnet-descriptor accelnet-setup-descriptor
    EXPORT AccelNetDescriptorsTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/"
    DESTINATION include OPTIONAL)
install(EXPORT AccelNetDescriptorsTargets
    NAMESPACE AccelNetDescriptors::
    DESTINATION lib/cmake/AccelNetDescriptors)
configure_package_config_file(
    cmake/AccelNetDescriptorsConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetDescriptorsConfig.cmake"
    INSTALL_DESTINATION lib/cmake/AccelNetDescriptors)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetDescriptorsConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetDescriptorsConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/AccelNetDescriptorsConfigVersion.cmake"
    DESTINATION lib/cmake/AccelNetDescriptors)
