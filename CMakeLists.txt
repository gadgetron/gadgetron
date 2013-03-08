
cmake_minimum_required(VERSION 2.6)

# add additional FindDCMTK Module
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
find_package(DCMTK REQUIRED)

set(CMAKE_INSTALL_PREFIX $ENV{GADGETRON_HOME})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -D_REENTRANT -D_OSF_SOURCE")

set(GT_DICOM_LIBRARIES
    z
    ${DCMTK_dcmdata_LIBRARY}
    ${DCMTK_oflog_LIBRARY}
    ${DCMTK_ofstd_LIBRARY}
    m
    #rt
    pthread
    #nsl
)

# sanity check:
message(STATUS "DCMTK ${DCMTK_HOME}")
message(STATUS "Include: ${DCMTK_INCLUDE_DIRS}")
message(STATUS "Libs: ${GT_DICOM_LIBRARIES}")

INCLUDE_DIRECTORIES(
    ${ISMRMRD_INCLUDE_DIR}
    ${ACE_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${GADGETRON_INCLUDE_DIR}
    ${DCMTK_INCLUDE_DIRS}
)

set(FMRIF_GE_LIB gadgetronfmrifgelib)
ADD_LIBRARY(${FMRIF_GE_LIB} SHARED DicomFinishGadget.cpp DicomImageWriter.cpp)

TARGET_LINK_LIBRARIES(
    ${FMRIF_GE_LIB}
    ${ISMRMRD_LIBRARIES}
    ${ISMRMRD_XSD_LIBRARIES}
    optimized ${ACE_LIBRARIES}
    debug ${ACE_DEBUG_LIBRARY}
    ${GT_DICOM_LIBRARIES}
)

INSTALL(
    FILES DicomFinishGadget.h DicomImageWriter.h fmrifgelib_export.h
    DESTINATION include
)

INSTALL(TARGETS ${FMRIF_GE_LIB} DESTINATION lib)

INSTALL(FILES dicom.xml DESTINATION config)
