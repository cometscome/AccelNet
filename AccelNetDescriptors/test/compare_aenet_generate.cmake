if(NOT DEFINED AENET_GENERATE OR NOT DEFINED AENET_TRNSET2ASCII OR
   NOT DEFINED ACCELNET_GENERATE OR NOT DEFINED COMPARATOR OR
   NOT DEFINED CASE_DIR OR NOT DEFINED STRUCTURE OR NOT DEFINED WORK)
    message(FATAL_ERROR "aenet generate comparison variables are missing")
endif()

file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(COPY
    "${CASE_DIR}/generate.in"
    "${CASE_DIR}/Ti.fingerprint.stp"
    "${CASE_DIR}/O.fingerprint.stp"
    DESTINATION "${WORK}")
file(COPY_FILE "${STRUCTURE}" "${WORK}/structure.xsf")

execute_process(
    COMMAND "${AENET_GENERATE}" generate.in
    WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE AENET_RESULT
    OUTPUT_FILE "${WORK}/aenet-generate.log"
    ERROR_FILE "${WORK}/aenet-generate.err")
if(NOT AENET_RESULT EQUAL 0)
    message(FATAL_ERROR "aenet generate.x failed: ${AENET_RESULT}")
endif()

execute_process(
    COMMAND "${AENET_TRNSET2ASCII}" --raw aenet.train aenet.txt
    WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE ASCII_RESULT
    OUTPUT_FILE "${WORK}/trnset2ascii.log"
    ERROR_FILE "${WORK}/trnset2ascii.err")
if(NOT ASCII_RESULT EQUAL 0)
    message(FATAL_ERROR "aenet trnset2ASCII.x failed: ${ASCII_RESULT}")
endif()

execute_process(
    COMMAND "${ACCELNET_GENERATE}" accelnet.hex 2
        Ti.fingerprint.stp O.fingerprint.stp structure.xsf
    WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE ACCELNET_RESULT
    OUTPUT_FILE "${WORK}/accelnet-generate.log"
    ERROR_FILE "${WORK}/accelnet-generate.err")
if(NOT ACCELNET_RESULT EQUAL 0)
    message(FATAL_ERROR "AccelNet descriptor generation failed: ${ACCELNET_RESULT}")
endif()

execute_process(
    COMMAND "${COMPARATOR}" aenet.txt accelnet.hex
    WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE COMPARE_RESULT
    OUTPUT_VARIABLE COMPARE_OUTPUT
    ERROR_VARIABLE COMPARE_ERROR)
if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR
        "aenet/AccelNet descriptor comparison failed: ${COMPARE_RESULT}\n"
        "${COMPARE_OUTPUT}${COMPARE_ERROR}")
endif()
message(STATUS "${COMPARE_OUTPUT}")
