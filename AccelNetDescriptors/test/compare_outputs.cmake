if(NOT DEFINED NEW_EXE OR NOT DEFINED ORIGINAL_EXE OR NOT DEFINED COMPARATOR OR
   NOT DEFINED INPUTS OR NOT DEFINED WORK)
    message(FATAL_ERROR "comparison variables are missing")
endif()
string(REPLACE "," ";" INPUT_LIST "${INPUTS}")
file(MAKE_DIRECTORY "${WORK}")
set(NEW_OUTPUT "${WORK}/new.hex")
set(ORIGINAL_OUTPUT "${WORK}/original.hex")
execute_process(
    COMMAND "${NEW_EXE}" "${NEW_OUTPUT}" ${INPUT_LIST}
    RESULT_VARIABLE NEW_RESULT)
if(NOT NEW_RESULT EQUAL 0)
    message(FATAL_ERROR "new descriptor executable failed: ${NEW_RESULT}")
endif()
execute_process(
    COMMAND "${ORIGINAL_EXE}" "${ORIGINAL_OUTPUT}" ${INPUT_LIST}
    RESULT_VARIABLE ORIGINAL_RESULT)
if(NOT ORIGINAL_RESULT EQUAL 0)
    message(FATAL_ERROR "original AccelNet descriptor executable failed: ${ORIGINAL_RESULT}")
endif()
execute_process(
    COMMAND "${COMPARATOR}" "${NEW_OUTPUT}" "${ORIGINAL_OUTPUT}"
    RESULT_VARIABLE COMPARE_RESULT)
if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR
        "descriptor files differ:\n  ${NEW_OUTPUT}\n  ${ORIGINAL_OUTPUT}")
endif()
message(STATUS "New package and untouched AccelNet descriptor files agree within strict tolerance")
