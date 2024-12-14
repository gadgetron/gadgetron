if (MKL_FOUND)
    add_definitions(-DHAVE_MKL)
    # These are needed to get the linking to work properly when
    # MKL is installed, but Armadillo is NOT using it.
    list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
    list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif ()

set( pingvin_cmr_header_files
    CmrCartesianKSpaceBinningCineGadget.h
    CmrParametricMappingGadget.h
    CmrParametricT1SRMappingGadget.h
    CmrParametricT2MappingGadget.h
    PureCmrCartesianKSpaceBinningCineGadget.h)

set( pingvin_cmr_src_files
    CmrCartesianKSpaceBinningCineGadget.cpp
    CmrParametricMappingGadget.cpp
    CmrParametricT1SRMappingGadget.cpp
    CmrParametricT2MappingGadget.cpp
    PureCmrCartesianKSpaceBinningCineGadget.cpp)

set( config_BinningCine_files
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning.xml
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_MultiSeries.xml)

source_group(config/BinningCine FILES ${config_BinningCine_files})

set( config_Mapping_files
    config/Mapping/CMR_2DT_T1Mapping_SASHA.xml)

source_group(config/Mapping FILES ${config_Mapping_files})

set( config_cmr_files
    ${config_BinningCine_files}
    ${config_Mapping_files})

if (BUILD_PYTHON_SUPPORT)

    message("Build python cmr gadgets ... ")

    include_directories(${CMAKE_SOURCE_DIR}/toolboxes/python)

    set( pingvin_cmr_header_files ${pingvin_cmr_header_files} CmrRealTimeLAXCineAIAnalysisGadget.h)
    set( pingvin_cmr_src_files ${pingvin_cmr_src_files} CmrRealTimeLAXCineAIAnalysisGadget.cpp)

    set( config_LandmarkDetection_files
        config/LandmarkDetection/CMR_RTCine_LAX_AI.xml
        config/LandmarkDetection/CMR_RTCine_Recon.xml
        config/LandmarkDetection/CMR_Image_Chain_RTCine_LAX_AI.xml
        config/LandmarkDetection/stream_image_array.xml
        )

    set( python_LandmarkDetection_files
        config/LandmarkDetection/gadgetron_cmr_landmark_detection.py
        config/LandmarkDetection/gadgetron_cmr_landmark_detection_util.py
        )

    source_group(config/LandmarkDetection FILES ${config_LandmarkDetection_files} ${python_LandmarkDetection_files})

    set( config_cmr_files
        ${config_cmr_files}
        ${config_LandmarkDetection_files}
        ${python_LandmarkDetection_files} )

    install(FILES ${python_LandmarkDetection_files} DESTINATION ${PINGVIN_INSTALL_PYTHON_MODULE_PATH}/cmr_lax_landmark_detection COMPONENT main)
endif ()

add_library(pingvin_cmr SHARED
    ${pingvin_cmr_header_files}
    ${pingvin_cmr_src_files}
    ${config_cmr_files} )

set_target_properties(pingvin_cmr PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_cmr
    pingvin_core
    pingvin_mricore
    pingvin_toolbox_log
    pingvin_toolbox_cpucore
    pingvin_toolbox_cpufft
    pingvin_toolbox_cpuklt
    pingvin_toolbox_mri_core
    pingvin_toolbox_cpuoperator
    pingvin_toolbox_image_analyze_io
    pingvin_toolbox_cmr
)

if (BUILD_PYTHON_SUPPORT)
    target_link_libraries(pingvin_cmr pingvin_toolbox_python)
endif ()

install(FILES ${pingvin_cmr_header_files} DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${config_cmr_files} DESTINATION ${PINGVIN_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS pingvin_cmr
        EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)