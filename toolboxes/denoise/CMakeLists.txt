add_library(pingvin_toolbox_denoise SHARED
        non_local_means.cpp
        non_local_means.h
        non_local_bayes.cpp non_local_bayes.h)

set_target_properties(pingvin_toolbox_denoise PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_denoise
        pingvin_toolbox_log
        pingvin_toolbox_cpucore
        pingvin_toolbox_cpucore_math
        armadillo
        )
target_include_directories(pingvin_toolbox_denoise
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

install(TARGETS pingvin_toolbox_denoise
		EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        non_local_means.h
        non_local_bayes.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)