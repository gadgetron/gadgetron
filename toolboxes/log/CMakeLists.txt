add_library(pingvin_toolbox_log SHARED log.cpp)
target_include_directories(pingvin_toolbox_log
		PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
set_target_properties(pingvin_toolbox_log PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

install(TARGETS pingvin_toolbox_log
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES log.h DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)