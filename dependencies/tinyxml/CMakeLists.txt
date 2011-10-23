cmake_minimum_required(VERSION 2.6)

#set(CMAKE_OSX_ARCHITECTURES "i386")

project(TinyXML)

#set(CMAKE_CXX_FLAGS "-m32 -Wall")
set(CMAKE_CXX_FLAGS "-Wall")

add_library(tinyxml SHARED tinyxmlerror.cpp tinyxmlparser.cpp tinyxml.cpp)

install(TARGETS tinyxml DESTINATION lib)


