# DCMTK-necessary preprocessor flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -D_REENTRANT -D_OSF_SOURCE")

find_package(Boost COMPONENTS date_time REQUIRED)

if (WIN32)
    set(GT_DICOM_LIBRARIES ${DCMTK_dcmdata_LIBRARY} ${DCMTK_oflog_LIBRARY} ${DCMTK_ofstd_LIBRARY})
else ()
    set(GT_DICOM_LIBRARIES
        ${DCMTK_dcmdata_LIBRARY}
        ${DCMTK_oflog_LIBRARY}
        ${DCMTK_ofstd_LIBRARY}
        z m pthread)
    if (APPLE)
        list(APPEND GT_DICOM_LIBRARIES iconv)
    endif ()
endif ()



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

target_link_libraries(
    gadgetron_dicom
    gadgetron_core
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore 
    ${GT_DICOM_LIBRARIES} )

install(
    FILES ${gadgetron_dicom_header_files} 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_dicom DESTINATION lib COMPONENT main)

install(FILES ${gadgetron_dicom_config_files}  DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
