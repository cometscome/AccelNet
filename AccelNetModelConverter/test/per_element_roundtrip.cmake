file(REMOVE_RECURSE "${OUTPUT}")
execute_process(COMMAND "${EXE}" n2p2-to-accelnet "${FIXTURE}" "${OUTPUT}/accelnet-1"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Per-element n2p2-to-AccelNet conversion failed")
endif()
execute_process(COMMAND "${EXE}" accelnet-to-n2p2 "${OUTPUT}/n2p2"
    "${OUTPUT}/accelnet-1/H.nn.ascii" "${OUTPUT}/accelnet-1/O.nn.ascii"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Per-element AccelNet-to-n2p2 conversion failed")
endif()
execute_process(COMMAND "${VALIDATOR}" "${OUTPUT}/n2p2" "${DEPTH_FIXTURE}"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Converted per-element n2p2 model failed numerical validation")
endif()
execute_process(COMMAND "${EXE}" n2p2-to-accelnet "${OUTPUT}/n2p2" "${OUTPUT}/accelnet-2"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Per-element round-trip reload failed")
endif()
foreach(element H O)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${OUTPUT}/accelnet-1/${element}.nn.ascii"
        "${OUTPUT}/accelnet-2/${element}.nn.ascii" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Per-element round-trip network differs for ${element}")
    endif()
endforeach()
