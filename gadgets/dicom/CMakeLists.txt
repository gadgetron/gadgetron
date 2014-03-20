# DCMTK-necessary preprocessor flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CONFIG_H -D_REENTRANT -D_OSF_SOURCE")

if (WIN32)
    set(GT_DICOM_LIBRARIES
        #z
        ${DCMTK_dcmdata_LIBRARY}
        ${DCMTK_oflog_LIBRARY}
        ${DCMTK_ofstd_LIBRARY}
        #m
        #rt
        #nsl
        #pthread
        )
else (WIN32)
    set(GT_DICOM_LIBRARIES
        z
        ${DCMTK_dcmdata_LIBRARY}
        ${DCMTK_oflog_LIBRARY}
        ${DCMTK_ofstd_LIBRARY}
        m
        #rt
        #nsl
        pthread)
endif (WIN32)

# sanity check:
#message("DCMTK ${DCMTK_HOME}")
#message("Include: ${DCMTK_INCLUDE_DIRS}")
#message("Libraries: ${GT_DICOM_LIBRARIES}")

include_directories(
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core    # for GadgetIsmrmrdReadWrite.h
    ${DCMTK_INCLUDE_DIRS})

set(GT_DICOM_LIB gadgetron_dicom)
add_library(${GT_DICOM_LIB} SHARED
    DicomFinishGadget.cpp
    DicomImageWriter.cpp
    ${ISMRMRD_XSD_SOURCE})

target_link_libraries(
    ${GT_DICOM_LIB}
    cpucore 
    ${ISMRMRD_LIBRARIES}
    optimized ${ACE_LIBRARIES}
    debug ${ACE_DEBUG_LIBRARY}
    ${XERCESC_LIBRARIES}
    ${GT_DICOM_LIBRARIES})

install(
    FILES DicomFinishGadget.h DicomImageWriter.h gadgetron_dicom_export.h
    DESTINATION include)

install(TARGETS ${GT_DICOM_LIB} DESTINATION lib)

install(FILES dicom.xml DESTINATION config)
