if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CMR__)
endif ()

if (MKL_FOUND)
    add_definitions(-DHAVE_MKL)
    # These are needed to get the linking to work properly when
    # MKL is installed, but Armadillo is NOT using it.
    list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
    list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif ()



set( gadgetron_cmr_header_files 
    gadgetron_cmr_export.h 
    CmrCartesianKSpaceBinningCineGadget.h 
    CmrParametricMappingGadget.h 
    CmrParametricT1SRMappingGadget.h 
    CmrParametricT2MappingGadget.h
    PureCmrCartesianKSpaceBinningCineGadget.h)

set( gadgetron_cmr_src_files 
    CmrCartesianKSpaceBinningCineGadget.cpp 
    CmrParametricMappingGadget.cpp 
    CmrParametricT1SRMappingGadget.cpp 
    CmrParametricT2MappingGadget.cpp
    PureCmrCartesianKSpaceBinningCineGadget.cpp)

set( config_BinningCine_files
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning.xml
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_Cloud.xml 
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_MultiSeries.xml
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_MultiSeries_Cloud.xml )

source_group(config/BinningCine FILES ${config_BinningCine_files})

set( config_Mapping_files 
    config/Mapping/CMR_2DT_T1Mapping_SASHA.xml)

source_group(config/Mapping FILES ${config_Mapping_files})

set( config_cmr_files 
    ${config_BinningCine_files} 
    ${config_Mapping_files})

if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)

    message("Build python cmr gadgets ... ")

    set( gadgetron_cmr_header_files ${gadgetron_cmr_header_files} CmrRealTimeLAXCineAIAnalysisGadget.h)
    set( gadgetron_cmr_src_files ${gadgetron_mricore_src_files} CmrRealTimeLAXCineAIAnalysisGadget.cpp)
    set( config_cmr_files ${config_cmr_files} config/LandmarkDetection/CMR_RTCine_LAX_AI.xml)

    set( python_LandmarkDetection_files 
        config/LandmarkDetection/gadgetron_cmr_landmark_detection.py
        config/LandmarkDetection/gadgetron_cmr_landmark_detection_util.py
        )

    install(FILES ${python_LandmarkDetection_files} DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)
endif ()

add_library(gadgetron_cmr SHARED 
    ${gadgetron_cmr_header_files} 
    ${gadgetron_cmr_src_files})

set_target_properties(gadgetron_cmr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_cmr 
    gadgetron_core
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuklt
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_image_analyze_io
    gadgetron_toolbox_cmr
)

install(FILES ${gadgetron_cmr_header_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${config_cmr_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_cmr
        EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)