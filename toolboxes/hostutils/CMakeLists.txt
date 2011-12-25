if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_HOSTUTILS__)
endif (WIN32)


include_directories(${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_library(hostutils SHARED parameterparser.cpp)

install(TARGETS hostutils DESTINATION lib)
install(FILES parameterparser.h hostutils_export.h FileInfo.h DESTINATION include)
