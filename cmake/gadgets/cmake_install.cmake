# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/gadgets

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/examples/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/interventional_mri/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cmr/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/epi/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/fatwater/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/T1/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/pmri/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/radial/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/hyper/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/gpu/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_noncartesian/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/spiral/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/grappa/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/python/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/dicom/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/cartesian/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/moco/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/plplot/cmake_install.cmake")

endif()

