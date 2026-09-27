# CMake integration for the LAMMPS ACCELNET package.

set(ACCELNET_DIR "" CACHE PATH "Path to the AccelNet source tree")

find_library(ACCELNET_LIBRARY NAMES accelnet
  HINTS "${ACCELNET_DIR}/build/lib" "${ACCELNET_DIR}/lib")
find_library(ACCELNET_DESCRIPTORS_LIBRARY NAMES AccelNetDescriptors
  HINTS "${ACCELNET_DIR}/build/lib" "${ACCELNET_DIR}/lib")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACCELNET DEFAULT_MSG
  ACCELNET_LIBRARY ACCELNET_DESCRIPTORS_LIBRARY)

if(NOT ACCELNET_FOUND)
  message(FATAL_ERROR "AccelNet was not found; set ACCELNET_DIR to the AccelNet source tree")
endif()

enable_language(Fortran)

target_link_directories(lammps PUBLIC ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
target_link_libraries(lammps PRIVATE
  "${ACCELNET_LIBRARY}"
  "${ACCELNET_DESCRIPTORS_LIBRARY}"
)
if(NOT PKG_GPU)
  target_link_libraries(lammps PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()

mark_as_advanced(ACCELNET_LIBRARY ACCELNET_DESCRIPTORS_LIBRARY)

# The optional /gpu style uses CUDA lib/gpu neighbors and Fortran OpenMP target.
if(PKG_GPU)
  set(ACCELNET_TARGET_FLAGS "-mp=gpu" CACHE STRING "NVHPC flags for final GPU link")
  find_library(ACCELNET_TARGET_LIBRARY NAMES accelnet_target
    HINTS "${ACCELNET_DIR}/build/lib" "${ACCELNET_DIR}/lib" REQUIRED)
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC" OR NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    message(FATAL_ERROR "ACCELNET GPU currently requires NVHPC C++ and Fortran compilers")
  endif()
  if(NOT GPU_API STREQUAL "cuda" AND NOT GPU_API STREQUAL "CUDA")
    message(FATAL_ERROR "ACCELNET GPU requires GPU_API=cuda")
  endif()
  if(NOT GPU_PREC STREQUAL "double" AND NOT GPU_PREC STREQUAL "DOUBLE")
    message(FATAL_ERROR "ACCELNET GPU requires GPU_PREC=double")
  endif()
  separate_arguments(_accelnet_target_flags NATIVE_COMMAND "${ACCELNET_TARGET_FLAGS}")
  target_link_libraries(lammps PRIVATE "${ACCELNET_TARGET_LIBRARY}" "${ACCELNET_LIBRARY}" "${ACCELNET_DESCRIPTORS_LIBRARY}")
  target_link_options(lammps PUBLIC ${_accelnet_target_flags} -fortranlibs)
  # LAMMPS appends this list again later. Let the NVHPC driver order its own
  # runtimes; placing libnvomp before the offload runtime makes device count 0.
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
endif()
