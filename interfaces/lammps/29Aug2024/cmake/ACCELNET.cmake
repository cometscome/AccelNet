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
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

mark_as_advanced(ACCELNET_LIBRARY ACCELNET_DESCRIPTORS_LIBRARY)
