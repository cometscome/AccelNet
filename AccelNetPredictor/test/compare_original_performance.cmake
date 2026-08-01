execute_process(
    COMMAND "${NEW_EXE}" 2 "${TI_SETUP}" "${O_SETUP}" "${TI_ASCII}" "${O_ASCII}" "${STRUCTURE}" "${REPEATS}"
    RESULT_VARIABLE new_result OUTPUT_VARIABLE new_output ERROR_VARIABLE new_error)
if(NOT new_result EQUAL 0)
    message(FATAL_ERROR "new predictor benchmark failed: ${new_error}")
endif()
execute_process(
    COMMAND "${ORIGINAL_EXE}" "${TI_BINARY}" "${O_BINARY}" "${STRUCTURE}" "${REPEATS}"
    RESULT_VARIABLE original_result OUTPUT_VARIABLE original_output ERROR_VARIABLE original_error)
if(NOT original_result EQUAL 0)
    message(FATAL_ERROR "original AccelNet benchmark failed: ${original_error}")
endif()
string(REGEX MATCH "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE[ \t]+([^ \n\r]+)" ignored "${new_output}")
set(new_time "${CMAKE_MATCH_1}")
string(REGEX MATCH "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE[ \t]+([^ \n\r]+)" ignored "${original_output}")
set(original_time "${CMAKE_MATCH_1}")
if(new_time STREQUAL "" OR original_time STREQUAL "")
    message(FATAL_ERROR "could not parse benchmark output")
endif()
message(STATUS "New predictor energy+force API: ${new_time} s/structure")
message(STATUS "Original AccelNet API:         ${original_time} s/structure")
if(new_time GREATER original_time)
    message(FATAL_ERROR "new structure API is slower than original AccelNet")
endif()
