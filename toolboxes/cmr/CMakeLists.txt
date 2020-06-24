if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CMR__)
endif ()

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${ARMADILLO_INCLUDE_DIRS} 
    ${ACE_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_BINARY_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/ffd
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
  )

set(cmr_header_fiels cmr_export.h 
                    cmr_kspace_binning.h
                    cmr_time_stamp.h 
                    cmr_motion_correction.h 
                    cmr_parametric_mapping.h 
                    cmr_t1_mapping.h 
                    cmr_t2_mapping.h 
                    cmr_spirit_recon.h 
                    cmr_strain_analysis.h 
                    cmr_radial_thickening.h 
                    cmr_analytical_strain.h 
                    cmr_file_and_directory_handling.h 
                    cmr_image_container_util.h 
                    cmr_ismrmrd_util.h 
                    )

set(cmr_src_fiels cmr_kspace_binning.cpp 
                cmr_time_stamp.cpp 
                cmr_motion_correction.cpp 
                cmr_parametric_mapping.cpp 
                cmr_t1_mapping.cpp 
                cmr_t2_mapping.cpp 
                cmr_spirit_recon.cpp 
                cmr_strain_analysis.cpp 
                cmr_radial_thickening.cpp 
                cmr_analytical_strain.cpp 
                cmr_file_and_directory_handling.cpp 
                cmr_image_container_util.cpp 
                cmr_ismrmrd_util.cpp 
                )

add_library(gadgetron_toolbox_cmr SHARED 
            ${cmr_header_fiels} 
            ${cmr_src_fiels} )

set_target_properties(gadgetron_toolbox_cmr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cmr
                    ${ISMRMRD_LIBRARIES}
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_image_analyze_io 
                    gadgetron_toolbox_mri_core 
                    gadgetron_toolbox_cpudwt 
                    gadgetron_toolbox_cpuoperator 
                    gadgetron_toolbox_cpu_image )

install(TARGETS gadgetron_toolbox_cmr DESTINATION lib COMPONENT main)

install(FILES ${cmr_header_fiels} 
      DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
