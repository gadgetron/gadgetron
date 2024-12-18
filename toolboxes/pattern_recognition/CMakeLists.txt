set(pr_header_files pr_kmeans.h )

set(pr_src_files pr_kmeans.cpp )

add_library(pingvin_toolbox_pr SHARED
            ${pr_header_files}
            ${pr_src_files} )

set_target_properties(pingvin_toolbox_pr PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_pr
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_log
                    pingvin_toolbox_cpucore_math
                    pingvin_toolbox_cpufft
                    pingvin_toolbox_cpuklt
                    pingvin_toolbox_image_analyze_io
                    pingvin_toolbox_cpudwt
                    pingvin_toolbox_cpu_image )
target_include_directories(pingvin_toolbox_pr
        PUBLIC
        $<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

install(TARGETS pingvin_toolbox_pr
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES ${pr_header_files}
      DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)