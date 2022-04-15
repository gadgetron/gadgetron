# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/libgadgetron_toolbox_cpucore.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/libgadgetron_toolbox_cpucore.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/libgadgetron_toolbox_cpucore.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpucore.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/../NDArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/../complext.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/../GadgetronException.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/../GadgetronTimer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/cpucore_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray_converter.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray_iterators.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray_utils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDArray_fileio.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho2DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho2DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho3DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho3DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho4DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho4DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho5DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho5DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho6DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho6DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho7DArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/ho7DArray.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoMatrix.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoMatrix.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDPoint.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDBoundaryHandler.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDBoundaryHandler.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDInterpolator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDInterpolatorNearestNeighbor.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDInterpolatorLinear.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/hoNDInterpolatorBSpline.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/image/hoNDImage.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/image/hoNDImage.hxx"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/image/hoNDImageContainer2D.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/algorithm/hoNDBSpline.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/cpu/algorithm/hoNDBSpline.hxx"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils/cmake_install.cmake")

endif()

