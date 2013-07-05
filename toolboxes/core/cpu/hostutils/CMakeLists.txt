if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_HOSTUTILS__)
endif (WIN32)

add_library(hostutils SHARED 
  parameterparser.cpp
  )

install(TARGETS 
  hostutils 
  DESTINATION lib)

install(FILES 
  hostutils_export.h 
  parameterparser.h 
  url_encode.h 
  FileInfo.h 
  DESTINATION include)
