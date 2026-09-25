cmake_minimum_required(VERSION 3.20)
project(AccelNetModelConverterFortran VERSION 1.0.1 LANGUAGES Fortran)

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
