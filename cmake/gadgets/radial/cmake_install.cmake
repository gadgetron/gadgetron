# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/gadgets/radial

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/RadialPhaseCorrectionGadget.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/gadgetron_radial_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/gpuRadialPrepGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/gpuRadialSensePrepGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/gpuRadialSpiritPrepGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/gpuRetroGatedSensePrepGadget.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial/libgadgetron_gpuradial.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial/libgadgetron_gpuradial.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/usr/local/lib:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial/libgadgetron_gpuradial.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/pmri:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes:/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/usr/local/lib:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_gpuradial.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial/config/cmake_install.cmake")

endif()

