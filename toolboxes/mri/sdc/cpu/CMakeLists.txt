if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUSDC__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

add_library(gadgetron_toolbox_cpusdc SHARED
  hoSDC.h
  hoSDC.cpp
  hoSDC_conv.h
  hoSDC_conv.cpp
  hoSDC_kernel.h
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
)

install(TARGETS gadgetron_toolbox_cpusdc
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
  hoSDC.h
  hoSDC_conv.h
  hoSDC_kernel.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH}
  COMPONENT main
)