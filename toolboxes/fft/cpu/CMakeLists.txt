add_library(pingvin_toolbox_cpufft SHARED
        hoNDFFT.h
        hoNDFFT.cpp
        )

set_target_properties(pingvin_toolbox_cpufft PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_include_directories(
        pingvin_toolbox_cpufft
        PUBLIC
        $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(pingvin_toolbox_cpufft
        pingvin_toolbox_cpucore
        pingvin_toolbox_log
        pingvin_toolbox_cpucore_math
        FFTW
        armadillo
        )



install(TARGETS pingvin_toolbox_cpufft
		EXPORT pingvin-export
		LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        hoNDFFT.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)