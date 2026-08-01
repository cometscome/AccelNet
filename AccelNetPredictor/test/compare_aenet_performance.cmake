execute_process(
    COMMAND "${NEW_EXE}" 2 "${TI_SETUP}" "${O_SETUP}" "${TI_NETWORK}" "${O_NETWORK}" "${STRUCTURE}" "${REPEATS}"
    RESULT_VARIABLE new_result OUTPUT_VARIABLE new_output ERROR_VARIABLE new_error)
if(NOT new_result EQUAL 0)
    message(FATAL_ERROR "new predictor benchmark failed: ${new_error}")
endif()
execute_process(
    COMMAND "${AENET_EXE}" "${TI_NETWORK}" "${O_NETWORK}" "${STRUCTURE}" "${REPEATS}"
    RESULT_VARIABLE aenet_result OUTPUT_VARIABLE aenet_output ERROR_VARIABLE aenet_error)
if(NOT aenet_result EQUAL 0)
    message(FATAL_ERROR "ænet benchmark failed: ${aenet_error}")
endif()

function(read_timing output label result)
    string(REGEX MATCH "${label}[ \t]+([^ \n\r]+)" ignored "${output}")
    if(CMAKE_MATCH_1 STREQUAL "")
        message(FATAL_ERROR "could not parse ${label}")
    endif()
    set(${result} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()

read_timing("${new_output}" "ENERGY_API_SECONDS_PER_STRUCTURE" new_energy)
read_timing("${new_output}" "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE" new_force)
read_timing("${aenet_output}" "ENERGY_API_SECONDS_PER_STRUCTURE" aenet_energy)
read_timing("${aenet_output}" "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE" aenet_force)
message(STATUS "New / ænet energy API:       ${new_energy} / ${aenet_energy} s/structure")
message(STATUS "New / ænet energy+force API: ${new_force} / ${aenet_force} s/structure")
if(new_energy GREATER aenet_energy)
    message(FATAL_ERROR "new energy structure API is slower than ænet")
endif()
if(new_force GREATER aenet_force)
    message(FATAL_ERROR "new energy+force structure API is slower than ænet")
endif()
