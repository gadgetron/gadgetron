add_library(pingvin_toolbox_cpuklt SHARED
  hoNDKLT.h
  hoNDKLT.cpp
  )

set_target_properties(pingvin_toolbox_cpuklt PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_cpuklt
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_log
                    pingvin_toolbox_cpucore_math
                    )



target_include_directories(
        pingvin_toolbox_cpuklt
        PUBLIC
        $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>

)

install(TARGETS pingvin_toolbox_cpuklt
	EXPORT pingvin-export
		LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
  hoNDKLT.h
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)
