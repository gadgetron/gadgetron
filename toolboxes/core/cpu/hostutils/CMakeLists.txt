if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_HOSTUTILS__)
endif (WIN32)

add_library(gadgetron_toolbox_hostutils SHARED 
  parameterparser.cpp
  )

set_target_properties(gadgetron_toolbox_hostutils PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

install(TARGETS 
  gadgetron_toolbox_hostutils 
  DESTINATION lib COMPONENT toolboxes-lib)

install(FILES 
  hostutils_export.h 
  parameterparser.h 
  url_encode.h 
  FileInfo.h 
  DESTINATION include COMPONENT toolboxes-dev)
