# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/javeda2/anaconda3/envs/gadgetron")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/generalOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/linearOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/identityOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/diagonalOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/diagonalSumOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/encodingOperatorContainer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/multiplicationOperatorContainer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/FFTOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/imageOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/encodedImageOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/partialDerivativeOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/convolutionOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/laplaceOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/downsampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/upsampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/tvPicsOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/subsetOperator.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu/cmake_install.cmake")

endif()

