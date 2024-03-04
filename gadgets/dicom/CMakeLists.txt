# DCMTK-necessary preprocessor flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -D_REENTRANT -D_OSF_SOURCE")

# set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost COMPONENTS date_time REQUIRED)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()

if (WIN32)
    set(GT_DICOM_LIBRARIES cmr i2d ijg8 dcmfg)
else ()
    set(GT_DICOM_LIBRARIES
            dcmdata
            oflog
            ofstd
        z m pthread)
    if (APPLE)
        list(APPEND GT_DICOM_LIBRARIES iconv)
    endif ()
endif ()

# sanity check:
message("DCMTK ${DCMTK_HOME}")
message("Include: ${DCMTK_INCLUDE_DIRS}")
message("Libraries: ${GT_DICOM_LIBRARIES}")

include_directories(
        ${CMAKE_SOURCE_DIR}/gadgets/mri_core    # for GadgetIsmrmrdReadWrite.h
        ${DCMTK_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${CMAKE_SOURCE_DIR}/gadgets/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
        ${CMAKE_SOURCE_DIR}/apps/gadgetron
        ${CMAKE_SOURCE_DIR}/apps/matlab
        ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
        ${ARMADILLO_INCLUDE_DIRS}
        )

set(gadgetron_dicom_header_files 
    gadgetron_dicom_export.h 
    DicomFinishGadget.h 
    DicomImageWriter.h 
    dicom_ismrmrd_utility.h )

set(gadgetron_dicom_src_files 
    DicomFinishGadget.cpp 
    DicomImageWriter.cpp 
    dicom_ismrmrd_utility.cpp )

set(gadgetron_dicom_config_files dicom.xml )

add_library(gadgetron_dicom SHARED
            ${gadgetron_dicom_header_files}
            ${gadgetron_dicom_src_files} 
            ${gadgetron_dicom_config_files} )

set_target_properties(gadgetron_dicom PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_dicom 
    gadgetron_core
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    ISMRMRD::ISMRMRD
    ${GT_DICOM_LIBRARIES} )

if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_dicom gadgetron_toolbox_cpucore_math )
endif()

install(
    FILES ${gadgetron_dicom_header_files} 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_dicom
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES ${gadgetron_dicom_config_files}  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
