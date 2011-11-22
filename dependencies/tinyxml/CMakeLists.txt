cmake_minimum_required(VERSION 2.6)

if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_TINYXML__)
ADD_DEFINITIONS(-DTIXML_USE_STL)
endif (WIN32)

project(TinyXML)

add_library(tinyxml SHARED tinyxmlerror.cpp tinyxmlparser.cpp tinyxml.cpp)

install(TARGETS tinyxml DESTINATION lib)

install(FILES
tinyxml.h
gadgetron_tinyxml_export.h
DESTINATION include)
