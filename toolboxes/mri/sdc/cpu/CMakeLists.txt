add_library(pingvin_toolbox_cpusdc SHARED
  hoSDC.cpp
)

set_target_properties(pingvin_toolbox_cpusdc PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_include_directories(pingvin_toolbox_cpusdc
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(pingvin_toolbox_cpusdc
  pingvin_toolbox_cpucore
  pingvin_toolbox_log
  pingvin_toolbox_operator
  pingvin_toolbox_nfft
  pingvin_toolbox_cpunfft
)

install(TARGETS pingvin_toolbox_cpusdc
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH}
  COMPONENT main
)