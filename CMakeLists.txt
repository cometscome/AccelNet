cmake_minimum_required(VERSION 3.20)
project(AccelNet VERSION 1.0.1 LANGUAGES Fortran C)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

include(CTest)

option(BUILD_SHARED_LIBS "Build AccelNet libraries as shared libraries" OFF)
option(ACCELNET_BUILD_REFERENCE_TESTS
    "Build tests against an external original AccelNet source tree" OFF)

# Keep all release artifacts in predictable locations for a top-level build.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# The standalone projects historically locate these data trees as siblings.
# Set monorepo-aware defaults while still allowing users to override them.
set(ACCELNET_ORIGINAL_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNet-a76ec3a83c8b0246328e9cca4817cc0afc84c925"
    CACHE PATH "Path to the untouched original AccelNet source tree")
set(ACCELNET_PREDICTOR_GOLDEN_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../AccelNetTrainer.jl/test/data/fortran_predict"
    CACHE PATH "Existing immutable NN/XSF golden data")
set(TIO2_DATASET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../TiO2-xsf"
    CACHE PATH "Optional TiO2 XSF corpus for integration tests")

set(ACCELNET_DESCRIPTORS_BUILD_TESTING "${BUILD_TESTING}"
    CACHE BOOL "Build descriptor tests" FORCE)
set(ACCELNET_DESCRIPTORS_BUILD_REFERENCE
    "${ACCELNET_BUILD_REFERENCE_TESTS}"
    CACHE BOOL "Build descriptor reference tests" FORCE)
set(ACCELNET_PREDICTOR_BUILD_TESTING "${BUILD_TESTING}"
    CACHE BOOL "Build predictor tests" FORCE)

# Dependency order: ModelConverter -> Predictor -> Descriptors.
add_subdirectory(AccelNetDescriptors)
add_subdirectory(AccelNetPredictor)
add_subdirectory(AccelNetModelConverter)
