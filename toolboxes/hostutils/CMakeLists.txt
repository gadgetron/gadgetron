
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_library(hostutils SHARED parameterparser.cpp)

install(TARGETS hostutils DESTINATION lib)
install(FILES parameterparser.h DESTINATION include)
