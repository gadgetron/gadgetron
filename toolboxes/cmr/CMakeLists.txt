set(cmr_header_files
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
                    )

set(cmr_src_files cmr_kspace_binning.cpp
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
                )

add_library(pingvin_toolbox_cmr SHARED
            ${cmr_header_files}
            ${cmr_src_files} )

set_target_properties(pingvin_toolbox_cmr PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_cmr
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_log
                    pingvin_toolbox_cpucore_math
                    pingvin_toolbox_cpufft
                    pingvin_toolbox_cpuklt
                    pingvin_toolbox_image_analyze_io
                    pingvin_toolbox_mri_core
                    pingvin_toolbox_cpudwt
                    pingvin_toolbox_cpuoperator
                    pingvin_toolbox_cpu_image )

target_include_directories(pingvin_toolbox_cmr
        PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/toolboxes/ffd>)

install(TARGETS pingvin_toolbox_cmr
		EXPORT pingvin-export
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
		RUNTIME DESTINATION bin
		COMPONENT main
)

install(FILES ${cmr_header_files}
      DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)