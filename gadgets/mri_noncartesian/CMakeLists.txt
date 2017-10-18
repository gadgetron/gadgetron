if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRI_NONCARTESIAN__)
endif ()

find_package(HDF5 1.8 REQUIRED COMPONENTS C)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu    
		${CMAKE_SOURCE_DIR}/toolboxes/nfft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${HDF5_C_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
)   

set( gadgetron_mri_noncartesian_header_files
  GriddingReconGadget.h
	CPUGriddingReconGadget.h)

set( gadgetron_mri_noncartesian_src_files
  GriddingReconGadget.cpp
	CPUGriddingReconGadget.cpp) 

set( gadgetron_mri_noncartesian_config_files
	config/Generic_CPU_Gridding_Recon.xml
  config/Generic_Spiral.xml
  config/Generic_Spiral_SNR.xml)

add_library(gadgetron_mri_noncartesian SHARED 
  gadgetron_mri_noncartesian_export.h 
  ${gadgetron_mri_noncartesian_header_files} 
  ${gadgetron_mri_noncartesian_src_files}
  ${gadgetron_mri_noncartesian_config_files}
)

set_target_properties(gadgetron_mri_noncartesian PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mri_noncartesian 
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_mri_core
    gadgetron_mricore
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
		gadgetron_toolbox_cpunfft
    gadgetron_toolbox_mri_core
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_image_analyze_io
    gadgetron_toolbox_gpucore
    gadgetron_toolbox_gpusolvers
    gadgetron_toolbox_gpuoperators
    gadgetron_toolbox_gpuparallelmri
    ${ISMRMRD_LIBRARIES} 
    ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${OPTIMIZED_GADGET_LIBS}
    ${Boost_LIBRARIES}

)

install(FILES 
    gadgetron_mri_noncartesian_export.h
    ${gadgetron_mricore_header_files}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_mri_noncartesian DESTINATION lib COMPONENT main)

install(FILES ${gadgetron_mri_noncartesian_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
