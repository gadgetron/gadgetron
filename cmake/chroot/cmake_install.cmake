# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/chroot

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xscriptsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/chroot" TYPE PROGRAM FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/copy-cuda-lib.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/start-gadgetron.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/enter-chroot-env.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/gadgetron-dependency-query.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/siemens_to_ismrmrd.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/gadgetron_ismrmrd_client.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/gadgetron_ismrmrd_client_noise_summary.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/gt_alive.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/run-gadgetron-dependency-query.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/run-gadgetron_ismrmrd_client.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/run-gadgetron_ismrmrd_client_noise_summary.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/run-gt_alive.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/run-siemens_to_ismrmrd.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/start-env.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/start.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/mount_image.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/chroot/start-gadgetron-from-image.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/mount.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/stop.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/umount_image.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/install_chroot_image.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/clean_gadgetron_data.sh"
    "/home/javeda2/mrprogs/gadgetron_fork/chroot/nvidia-copy.sh"
    )
endif()

