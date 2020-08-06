if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUSDC__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

add_library(gadgetron_toolbox_cpusdc SHARED
  hoSDC.cpp
)

set_target_properties(gadgetron_toolbox_cpusdc PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(gadgetron_toolbox_cpusdc
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)
        
target_link_libraries(gadgetron_toolbox_cpusdc
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_hostutils
  gadgetron_toolbox_log
  gadgetron_toolbox_operator
  gadgetron_toolbox_nfft
  gadgetron_toolbox_cpunfft
)

install(TARGETS gadgetron_toolbox_cpusdc
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}
  COMPONENT main
)