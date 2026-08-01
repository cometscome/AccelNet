if(NOT DEFINED NEW_EXE OR NOT DEFINED LEGACY_EXE OR NOT DEFINED TI_SETUP OR
   NOT DEFINED O_SETUP OR NOT DEFINED INPUT OR NOT DEFINED WORK)
    message(FATAL_ERROR "setup CLI comparison variables are missing")
endif()
file(MAKE_DIRECTORY "${WORK}")
set(NEW_OUTPUT "${WORK}/setup.hex")
set(LEGACY_OUTPUT "${WORK}/legacy.hex")
execute_process(
    COMMAND "${NEW_EXE}" "${NEW_OUTPUT}" 2 "${TI_SETUP}" "${O_SETUP}" "${INPUT}"
    RESULT_VARIABLE NEW_RESULT)
if(NOT NEW_RESULT EQUAL 0)
    message(FATAL_ERROR "setup descriptor executable failed: ${NEW_RESULT}")
endif()
execute_process(
    COMMAND "${LEGACY_EXE}" "${LEGACY_OUTPUT}" "${INPUT}"
    RESULT_VARIABLE LEGACY_RESULT)
if(NOT LEGACY_RESULT EQUAL 0)
    message(FATAL_ERROR "legacy descriptor executable failed: ${LEGACY_RESULT}")
endif()
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${NEW_OUTPUT}" "${LEGACY_OUTPUT}"
    RESULT_VARIABLE COMPARE_RESULT)
if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR "setup-driven and direct descriptor files differ")
endif()
