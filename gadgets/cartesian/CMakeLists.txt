if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CARTESIAN__)
endif ()


add_library(gadgetron_cartesian SHARED 
  gadgetron_cartesian_export.h
  CartesianToGenericGadget.h
  CartesianToGenericGadget.cpp
)

set_target_properties(gadgetron_cartesian PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_cartesian gadgetron_toolbox_cpucore 
  gadgetron_core
  gadgetron_toolbox_log
  gadgetron_mricore
  gadgetron_toolbox_cpucore_math 
  )


install (FILES  gadgetron_cartesian_export.h
                CartesianToGenericGadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_cartesian
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)