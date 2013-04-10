find_package(Ismrmrd REQUIRED)

include_directories(
    ${MATLAB_INCLUDE_DIR}
    ${ACE_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
    ${CMAKE_SOURCE_DIR}/apps/gadgetron 
    ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
    ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
    ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
    ${CMAKE_SOURCE_DIR}/gadgets/core
)

if (UNIX)
    if (APPLE)
        SET(MATLAB_SUFFIX ".mexmaci64")
    else(APPLE)
        SET(MATLAB_SUFFIX ".mexglnxa64")
    endif(APPLE)
else(UNIX)
    SET(MATLAB_SUFFIX ".dll")
endif(UNIX)

add_library(gadgetronmatlab SHARED MatlabGadget.cpp)
target_link_libraries(
    gadgetronmatlab
    ${MATLAB_LIBRARIES}
    ${ISMRMRD_LIBRARIES}
    ${ISMRMRD_XSD_LIBRARIES}
    optimized ${ACE_LIBRARIES}
    debug ${ACE_DEBUG_LIBRARY}
)

#set(JAVA_MATLAB_SERVER_SRC "MatlabCommandServer.java")
#string(REPLACE "java" "class" JAVA_MATLAB_SERVER_CLASS ${JAVA_MATLAB_SERVER_SRC})
#add_custom_command(
#    OUTPUT ${JAVA_MATLAB_SERVER_CLASS}
#    COMMAND javac -cp "${MATLAB_JARS}" ${CMAKE_CURRENT_SOURCE_DIR}/${JAVA_MATLAB_SERVER_SRC}
#    COMMENT "Generating Matlab Command Server class" VERBATIM
#)
#add_custom_target(java_server ALL DEPENDS ${JAVA_MATLAB_SERVER_CLASS})

install(TARGETS gadgetronmatlab DESTINATION lib)
install(FILES BaseGadget.m scale.m accumulate_and_recon.m DESTINATION matlab)
install(FILES matlab.xml DESTINATION config)
#install(FILES ${JAVA_MATLAB_SERVER_CLASS} DESTINATION matlab)
