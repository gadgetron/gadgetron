cmake_minimum_required(VERSION 2.6)

if (WIN32)
ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_TINYXML__)
ADD_DEFINITIONS(-DTIXML_USE_STL)
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
endif (WIN32)

project(TinyXML)

add_library(tinyxml SHARED tinyxmlerror.cpp tinyxmlparser.cpp tinyxml.cpp)

install(TARGETS tinyxml DESTINATION lib)

install(FILES
tinyxml.h
gadgetron_tinyxml_export.h
DESTINATION include)
