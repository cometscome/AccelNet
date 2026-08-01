if(NOT DEFINED NEW_EXE OR NOT DEFINED ORIGINAL_EXE OR NOT DEFINED COMPARATOR OR NOT DEFINED VERSION OR
   NOT DEFINED TI_SETUP OR NOT DEFINED O_SETUP OR NOT DEFINED INPUT OR NOT DEFINED WORK)
    message(FATAL_ERROR "missing setup-version comparison argument")
endif()

file(MAKE_DIRECTORY "${WORK}")
set(new_output "${WORK}/new-v${VERSION}.hex")
set(original_output "${WORK}/original-v${VERSION}.hex")

execute_process(
    COMMAND "${NEW_EXE}" "${new_output}" 2 "${TI_SETUP}" "${O_SETUP}" "${INPUT}"
    RESULT_VARIABLE new_result)
if(NOT new_result EQUAL 0)
    message(FATAL_ERROR "setup-aware descriptor generator failed")
endif()

execute_process(
    COMMAND "${ORIGINAL_EXE}" "${original_output}" "--version=${VERSION}" "${INPUT}"
    RESULT_VARIABLE original_result)
if(NOT original_result EQUAL 0)
    message(FATAL_ERROR "original AccelNet descriptor generator failed")
endif()

execute_process(
    COMMAND "${COMPARATOR}" "${new_output}" "${original_output}"
    RESULT_VARIABLE comparison_result)
if(NOT comparison_result EQUAL 0)
    message(FATAL_ERROR "version ${VERSION} descriptor files exceed the comparison tolerance")
endif()
