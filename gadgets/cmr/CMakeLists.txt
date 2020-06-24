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

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
    ${CMAKE_SOURCE_DIR}/toolboxes/cmr
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/python
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_PATH}
    ${NUMPY_INCLUDE_DIRS}
)

set( gadgetron_cmr_header_files 
    gadgetron_cmr_export.h 
    CmrCartesianKSpaceBinningCineGadget.h 
    CmrParametricMappingGadget.h 
    CmrParametricT1SRMappingGadget.h 
    CmrParametricT2MappingGadget.h 
    CmrRealTimeLAXCineAIAnalysisGadget.h )

set( gadgetron_cmr_src_files 
    CmrCartesianKSpaceBinningCineGadget.cpp 
    CmrParametricMappingGadget.cpp 
    CmrParametricT1SRMappingGadget.cpp 
    CmrParametricT2MappingGadget.cpp 
    CmrRealTimeLAXCineAIAnalysisGadget.cpp )

set( config_BinningCine_files 
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning.xml
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_Cloud.xml 
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_MultiSeries.xml
    config/BinningCine/CMR_2DT_RTCine_KspaceBinning_MultiSeries_Cloud.xml 
    )

source_group(config/BinningCine FILES ${config_BinningCine_files})

set( config_Mapping_files 
    config/Mapping/CMR_2DT_T1Mapping_SASHA.xml
    )

source_group(config/Mapping FILES ${config_Mapping_files})

set( config_LandmarkDetection_files 
    config/LandmarkDetection/CMR_RTCine_LAX_AI.xml
    )

set( python_LandmarkDetection_files 
    config/LandmarkDetection/gadgetron_cmr_landmark_detection.py
    config/LandmarkDetection/gadgetron_cmr_landmark_detection_util.py
    )

source_group(config/LandmarkDetection FILES ${config_LandmarkDetection_files} ${python_LandmarkDetection_files})

set( config_cmr_files 
    ${config_BinningCine_files} 
    ${config_Mapping_files} 
    ${config_LandmarkDetection_files} 
    ${python_LandmarkDetection_files} )

add_library(gadgetron_cmr SHARED 
    ${gadgetron_cmr_header_files} 
    ${gadgetron_cmr_src_files}
    ${config_cmr_files}
)

set_target_properties(gadgetron_cmr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_cmr 
    gadgetron_gadgetbase
    gadgetron_mricore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuklt
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_image_analyze_io
    gadgetron_toolbox_cmr
    ${ISMRMRD_LIBRARIES} 
    ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${OPTIMIZED_GADGET_LIBS}
    ${Boost_LIBRARIES}
)

if(PYTHONLIBS_FOUND)
    target_link_libraries(gadgetron_cmr
                gadgetron_toolbox_python ${PYTHON_LIBRARIES} )
endif (PYTHONLIBS_FOUND)

install(FILES ${gadgetron_cmr_header_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${config_cmr_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
install(FILES ${python_LandmarkDetection_files} DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)

install(TARGETS gadgetron_cmr DESTINATION lib COMPONENT main)
