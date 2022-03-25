if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CMR__)
endif ()

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
                cmr_image_container_util.cpp 
                cmr_ismrmrd_util.cpp 
                )

add_library(gadgetron_toolbox_cmr SHARED 
            ${cmr_header_fiels} 
            ${cmr_src_fiels} )

set_target_properties(gadgetron_toolbox_cmr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cmr
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

target_include_directories(gadgetron_toolbox_cmr
        PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/ffd>)

install(TARGETS gadgetron_toolbox_cmr
		EXPORT gadgetron-export
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		RUNTIME DESTINATION bin
		COMPONENT main
)


install(FILES ${cmr_header_fiels} 
      DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)