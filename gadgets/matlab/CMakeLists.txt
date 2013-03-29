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

install(TARGETS gadgetronmatlab DESTINATION lib)

#install(FILES DESTINATION matlab)
install(FILES matlab.xml DESTINATION config)
