if(NOT DEFINED N2P2_SCALING OR NOT DEFINED ACCELNET_DUMP OR
   NOT DEFINED COMPARATOR OR NOT DEFINED CASE_DIR OR NOT DEFINED WORK OR
   NOT DEFINED CUTOFF_TYPE)
    message(FATAL_ERROR "missing n2p2 descriptor comparison argument")
endif()

file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(COPY "${CASE_DIR}/input.data" DESTINATION "${WORK}")
file(READ "${CASE_DIR}/input.nn.in" INPUT_NN)
string(REPLACE "@N2P2_CUTOFF_TYPE@" "${CUTOFF_TYPE}" INPUT_NN "${INPUT_NN}")
file(WRITE "${WORK}/input.nn" "${INPUT_NN}")

execute_process(
    COMMAND "${N2P2_SCALING}" 100
    WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE N2P2_RESULT
    OUTPUT_VARIABLE N2P2_OUTPUT
    ERROR_VARIABLE N2P2_ERROR)
if(NOT N2P2_RESULT EQUAL 0)
    message(FATAL_ERROR
        "nnp-scaling failed (${N2P2_RESULT})\nstdout:\n${N2P2_OUTPUT}\nstderr:\n${N2P2_ERROR}")
endif()
if(NOT EXISTS "${WORK}/function.data")
    message(FATAL_ERROR "nnp-scaling did not create function.data")
endif()

execute_process(
    COMMAND "${ACCELNET_DUMP}" "${WORK}" "${CASE_DIR}/structure.xsf"
            "${WORK}/accelnet.hex"
    RESULT_VARIABLE ACCELNET_RESULT
    OUTPUT_VARIABLE ACCELNET_OUTPUT
    ERROR_VARIABLE ACCELNET_ERROR)
if(NOT ACCELNET_RESULT EQUAL 0)
    message(FATAL_ERROR
        "AccelNet descriptor dump failed (${ACCELNET_RESULT})\nstdout:\n${ACCELNET_OUTPUT}\nstderr:\n${ACCELNET_ERROR}")
endif()

execute_process(
    COMMAND "${COMPARATOR}" "${WORK}/function.data" "${WORK}/accelnet.hex"
    RESULT_VARIABLE COMPARE_RESULT
    OUTPUT_VARIABLE COMPARE_OUTPUT
    ERROR_VARIABLE COMPARE_ERROR)
if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR
        "n2p2 descriptor comparison failed (${COMPARE_RESULT})\nstdout:\n${COMPARE_OUTPUT}\nstderr:\n${COMPARE_ERROR}")
endif()
message(STATUS "cutoff ${CUTOFF_TYPE}: ${COMPARE_OUTPUT}")
