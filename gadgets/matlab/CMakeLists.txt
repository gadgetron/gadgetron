cmake_minimum_required(VERSION 2.6)

# add additional FindMatlab Module
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

###############################################################
#Bootstrap search for libraries
# (We need to find cmake modules in Gadgetron)
###############################################################
find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
    $ENV{GADGETRON_HOME}/cmake
    /usr/local/gadgetron)

if (NOT GADGETRON_CMAKE_MODULES)
    MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found.
    Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

list(APPEND CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})
###############################################################

find_package(Ismrmrd REQUIRED)
find_package(Gadgetron REQUIRED)
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS thread system REQUIRED)
find_package(ACE REQUIRED)
find_package(Matlab REQUIRED)

set(CMAKE_INSTALL_PREFIX $ENV{GADGETRON_HOME})

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
