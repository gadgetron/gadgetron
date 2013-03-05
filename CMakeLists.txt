cmake_minimum_required(VERSION 2.6)

PROJECT(gt-tomography)

find_package(Gadgetron REQUIRED)
include_directories(  "${GADGETRON_HOME}/include")

LINK_DIRECTORIES(${GADGETRON_LIB_DIR})   
find_package(Boost REQUIRED)
find_package(CUDA 4.1)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
