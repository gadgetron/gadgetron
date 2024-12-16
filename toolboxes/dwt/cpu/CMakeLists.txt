add_library(pingvin_toolbox_cpudwt SHARED
    hoNDWavelet.h
    hoNDWavelet.cpp
    hoNDHarrWavelet.h
    hoNDHarrWavelet.cpp
    hoNDRedundantWavelet.h
    hoNDRedundantWavelet.cpp
    )
target_include_directories(pingvin_toolbox_cpudwt
        PUBLIC

        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_target_properties(pingvin_toolbox_cpudwt PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_cpudwt
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_log
                    pingvin_toolbox_cpucore_math
                    )

install(TARGETS pingvin_toolbox_cpudwt
        EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
    hoNDWavelet.h
    hoNDHarrWavelet.h
    hoNDRedundantWavelet.h
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)