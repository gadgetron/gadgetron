find_package(Ismrmrd REQUIRED)

include_directories(
  ${MATLAB_INCLUDE_DIR}
  #${GADGETRON_INCLUDE_DIR}
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

install(TARGETS gadgetronmatlab DESTINATION lib)
install(FILES MatlabGadget.h DESTINATION include)
install(FILES BaseGadget.m scale.m accumulate_and_recon.m mask_image.m DESTINATION matlab)
install(FILES matlab.xml DESTINATION config)
install(FILES ${JAVA_MATLAB_SERVER_CLASS} DESTINATION matlab)
