set(image_header_files
        morphology.h )

set(image_src_files morphology.cpp )

add_library(pingvin_toolbox_cpu_image SHARED
    ${image_header_files} ${image_src_files} )

set_target_properties(pingvin_toolbox_cpu_image PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_cpu_image
    pingvin_toolbox_log
    pingvin_toolbox_cpucore
    pingvin_toolbox_cpucore_math
    )

target_include_directories(pingvin_toolbox_cpu_image
        PUBLIC
        $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

install(TARGETS pingvin_toolbox_cpu_image
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
  ${image_header_files}
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)