cmake_minimum_required(VERSION 3.20)
project(AccelNetModelConverterFortran VERSION 1.1.1 LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(ACCELNET_PREDICTOR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetPredictor"
    CACHE PATH "AccelNetPredictor source directory")

if(TARGET AccelNetPredictor::AccelNetPredictor)
elseif(EXISTS "${ACCELNET_PREDICTOR_SOURCE_DIR}/CMakeLists.txt")
    set(ACCELNET_PREDICTOR_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    add_subdirectory("${ACCELNET_PREDICTOR_SOURCE_DIR}"
        "${CMAKE_CURRENT_BINARY_DIR}/AccelNetPredictor-build")
else()
    find_package(AccelNetPredictor CONFIG REQUIRED)
endif()

add_executable(accelnet-model-converter-fortran
    fortran/model_conversion.f90
    fortran/main.f90)
target_link_libraries(accelnet-model-converter-fortran PRIVATE
    AccelNetPredictor::AccelNetPredictor)
target_compile_options(accelnet-model-converter-fortran PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:-O2;-ffree-line-length-none>
    $<$<Fortran_COMPILER_ID:Intel,IntelLLVM>:-O2>)

install(TARGETS accelnet-model-converter-fortran RUNTIME DESTINATION bin)

include(CTest)
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/n2p2/input.nn")
    add_test(NAME fortran_converter_roundtrip
        COMMAND "${CMAKE_COMMAND}"
            -DEXE=$<TARGET_FILE:accelnet-model-converter-fortran>
            -DFIXTURE=${CMAKE_CURRENT_SOURCE_DIR}/test/fixtures/n2p2
            -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/fortran-roundtrip
            -P ${CMAKE_CURRENT_SOURCE_DIR}/test/fortran_roundtrip.cmake)
    if(TARGET test_n2p2_topology AND
       EXISTS "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/data/n2p2-per-element/input.nn")
        add_test(NAME fortran_converter_per_element_roundtrip
            COMMAND "${CMAKE_COMMAND}"
                -DEXE=$<TARGET_FILE:accelnet-model-converter-fortran>
                -DVALIDATOR=$<TARGET_FILE:test_n2p2_topology>
                -DFIXTURE=${ACCELNET_PREDICTOR_SOURCE_DIR}/test/data/n2p2-per-element
                -DDEPTH_FIXTURE=${ACCELNET_PREDICTOR_SOURCE_DIR}/test/data/n2p2-per-element-depth
                -DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/fortran-per-element-roundtrip
                -P ${CMAKE_CURRENT_SOURCE_DIR}/test/per_element_roundtrip.cmake)
        set_tests_properties(fortran_converter_per_element_roundtrip
            PROPERTIES LABELS "integration;n2p2-reference")
    endif()
endif()

# Synthetic checks are self-contained; an optional independent n2p2 executable
# adds an external oracle without making n2p2 a runtime/build dependency.
set(N2P2_REFERENCE_BENCHMARK "" CACHE FILEPATH "Independent serial n2p2-reference executable")
set(N2P2_PREDICT_EXAMPLES "" CACHE PATH "Optional n2p2/examples/nnp-predict directory")
if(BUILD_TESTING AND TARGET accelnet-n2p2-benchmark)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
        set(_n2p2_oracle)
        if(N2P2_REFERENCE_BENCHMARK)
            list(APPEND _n2p2_oracle --reference "${N2P2_REFERENCE_BENCHMARK}")
        endif()
        add_test(NAME n2p2_extended_multi_element
            COMMAND "${Python3_EXECUTABLE}" "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/check_n2p2_extended.py"
                ${_n2p2_oracle} --candidate $<TARGET_FILE:accelnet-n2p2-benchmark>
                --converter $<TARGET_FILE:accelnet-model-converter-fortran>
                --output "${CMAKE_CURRENT_BINARY_DIR}/extended-check")
        set_tests_properties(n2p2_extended_multi_element PROPERTIES TIMEOUT 300 LABELS "integration")
        if(N2P2_REFERENCE_BENCHMARK)
            add_test(NAME n2p2_high_g5_moments
                COMMAND "${Python3_EXECUTABLE}" "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/check_n2p2_high_moments.py"
                    ${_n2p2_oracle} --candidate $<TARGET_FILE:accelnet-n2p2-benchmark>
                    --converter $<TARGET_FILE:accelnet-model-converter-fortran>
                    --output "${CMAKE_CURRENT_BINARY_DIR}/high-moment-check")
            set_tests_properties(n2p2_high_g5_moments PROPERTIES TIMEOUT 300 LABELS "integration;n2p2-reference")
        endif()
        if(N2P2_REFERENCE_BENCHMARK AND N2P2_PREDICT_EXAMPLES)
            add_test(NAME n2p2_official_multi_element
                COMMAND "${Python3_EXECUTABLE}" "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/check_n2p2_official.py"
                    ${_n2p2_oracle} --candidate $<TARGET_FILE:accelnet-n2p2-benchmark>
                    --converter $<TARGET_FILE:accelnet-model-converter-fortran>
                    --examples "${N2P2_PREDICT_EXAMPLES}" --output "${CMAKE_CURRENT_BINARY_DIR}/official-check")
            set_tests_properties(n2p2_official_multi_element PROPERTIES TIMEOUT 900 LABELS "integration;n2p2-reference")
        endif()
    endif()
endif()

option(ACCELNET_TEST_N2P2_GPU "Run extended n2p2 checks on an actual GPU" OFF)
if(ACCELNET_TEST_N2P2_GPU AND BUILD_TESTING AND TARGET accelnet-n2p2-benchmark
   AND Python3_Interpreter_FOUND AND NOT ACCELNET_TARGET_SERIAL)
    if(N2P2_REFERENCE_BENCHMARK)
        add_test(NAME n2p2_high_g5_moments_gpu
            COMMAND "${Python3_EXECUTABLE}" "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/check_n2p2_high_moments.py"
                ${_n2p2_oracle} --candidate $<TARGET_FILE:accelnet-n2p2-benchmark>
                --converter $<TARGET_FILE:accelnet-model-converter-fortran> --backend gpu
                --output "${CMAKE_CURRENT_BINARY_DIR}/high-moment-gpu-check")
        set_tests_properties(n2p2_high_g5_moments_gpu PROPERTIES TIMEOUT 600
            LABELS "gpu;integration;n2p2-reference" ENVIRONMENT "OMP_TARGET_OFFLOAD=MANDATORY;OMP_NUM_THREADS=1")
    endif()
    add_test(NAME n2p2_extended_multi_element_gpu
        COMMAND "${Python3_EXECUTABLE}" "${ACCELNET_PREDICTOR_SOURCE_DIR}/test/check_n2p2_extended.py"
            ${_n2p2_oracle} --candidate $<TARGET_FILE:accelnet-n2p2-benchmark>
            --converter $<TARGET_FILE:accelnet-model-converter-fortran> --backend gpu
            --output "${CMAKE_CURRENT_BINARY_DIR}/extended-gpu-check")
    set_tests_properties(n2p2_extended_multi_element_gpu PROPERTIES TIMEOUT 600
        LABELS "gpu;integration" ENVIRONMENT "OMP_TARGET_OFFLOAD=MANDATORY;OMP_NUM_THREADS=1")
endif()
