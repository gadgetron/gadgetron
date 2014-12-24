find_package(Ismrmrd REQUIRED)

include_directories(${MATLAB_INCLUDE_DIR})

if (UNIX)
    if (APPLE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
        SET(MATLAB_SUFFIX ".mexmaci64")
    else(APPLE)
        SET(MATLAB_SUFFIX ".mexglnxa64")
    endif(APPLE)
else(UNIX)
    SET(MATLAB_SUFFIX ".dll")
endif(UNIX)

add_library(gadgetron_matlab SHARED gadgetron_matlab_export.h MatlabGadget.h MatlabGadget.cpp)
set_target_properties(gadgetron_matlab PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
    gadgetron_matlab
    gadgetron_gadgetbase
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    ${MATLAB_LIBRARIES}
    ${ISMRMRD_LIBRARIES}
    optimized ${ACE_LIBRARIES}
    debug ${ACE_DEBUG_LIBRARY}
)

if (UNIX)
    set(JAVA_MATLAB_SERVER_SRC "MatlabCommandServer.java")
    string(REPLACE "java" "class" JAVA_MATLAB_SERVER_CLASS ${JAVA_MATLAB_SERVER_SRC})
    set(JAVA_MATLAB_SERVER_CLASS "${CMAKE_CURRENT_BINARY_DIR}/${JAVA_MATLAB_SERVER_CLASS}")

    string(REPLACE ";" ":" MATLAB_UNIX_JARS "${MATLAB_JARS}")

    add_custom_command(
        OUTPUT ${JAVA_MATLAB_SERVER_CLASS}
        DEPENDS ${JAVA_MATLAB_SERVER_SRC}
        COMMAND javac -d ${CMAKE_CURRENT_BINARY_DIR} -cp "${MATLAB_UNIX_JARS}" ${CMAKE_CURRENT_SOURCE_DIR}/${JAVA_MATLAB_SERVER_SRC}
        COMMENT "Generating Matlab Command Server class" VERBATIM
    )
    add_custom_target(matlab_command_server ALL DEPENDS ${JAVA_MATLAB_SERVER_CLASS})
    install(FILES ${JAVA_MATLAB_SERVER_CLASS} DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
else(UNIX)
        MESSAGE( "Don't know how to build the Matlab Command Server class on Windows" )
endif(UNIX)

install(TARGETS gadgetron_matlab DESTINATION lib COMPONENT main)
install(FILES MatlabGadget.h gadgetron_matlab_export.h DESTINATION include COMPONENT main)
install(FILES BaseGadget.m scale.m accumulate_and_recon.m mask_image.m DESTINATION ${GADGETRON_INSTALL_MATLAB_PATH} COMPONENT main)
install(FILES matlab.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
