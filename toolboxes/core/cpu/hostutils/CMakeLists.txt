if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_HOSTUTILS__)
endif ()

add_library(gadgetron_toolbox_hostutils SHARED 
  parameterparser.cpp
  )

set_target_properties(gadgetron_toolbox_hostutils PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS 
  gadgetron_toolbox_hostutils 
  DESTINATION lib COMPONENT main)

install(FILES 
  hostutils_export.h 
  parameterparser.h 
  url_encode.h 
  FileInfo.h 
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
