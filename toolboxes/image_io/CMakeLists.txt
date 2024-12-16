set(image_io_header_files
        ImageIOBase.h
        ImageIOAnalyze.h)

set(image_io_src_files
        ImageIOBase.cpp
        ImageIOAnalyze.cpp)

add_library(pingvin_toolbox_image_analyze_io SHARED ${image_io_header_files} ${image_io_src_files})
set_target_properties(pingvin_toolbox_image_analyze_io PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})
target_link_libraries(pingvin_toolbox_image_analyze_io pingvin_toolbox_log pingvin_toolbox_cpucore Boost::filesystem)
target_include_directories(pingvin_toolbox_image_analyze_io
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)


install(TARGETS pingvin_toolbox_image_analyze_io
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        ${image_io_header_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)