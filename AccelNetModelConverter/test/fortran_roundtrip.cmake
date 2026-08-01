file(REMOVE_RECURSE "${OUTPUT}")
execute_process(COMMAND "${EXE}" n2p2-to-accelnet "${FIXTURE}" "${OUTPUT}/accelnet-1"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Fortran n2p2-to-AccelNet conversion failed")
endif()
execute_process(COMMAND "${EXE}" accelnet-to-n2p2 "${OUTPUT}/n2p2"
    "${OUTPUT}/accelnet-1/O.nn.ascii" "${OUTPUT}/accelnet-1/Ti.nn.ascii"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Fortran AccelNet-to-n2p2 conversion failed")
endif()
execute_process(COMMAND "${EXE}" n2p2-to-accelnet "${OUTPUT}/n2p2" "${OUTPUT}/accelnet-2"
    RESULT_VARIABLE result)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Fortran round-trip reload failed")
endif()
foreach(element O Ti)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${OUTPUT}/accelnet-1/${element}.nn.ascii"
        "${OUTPUT}/accelnet-2/${element}.nn.ascii" RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "Round-trip network differs for ${element}")
    endif()
endforeach()
