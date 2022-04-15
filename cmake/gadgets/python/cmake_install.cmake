# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/gadgets/python

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/python" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/accumulate_and_recon.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/bucket_recon.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/image_array_recon.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/image_array_recon_rtcine_plotting.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/passthrough.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/passthrough_array_image.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/pseudoreplicagather.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/remove_2x_oversampling.py"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/gadgets/rms_coil_combine.py"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/config" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/pseudoreplica.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/python_buckets.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/python_short.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/python_image_array_recon.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/python_passthrough.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/legacy/config/python_short.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/python/config/Generic_Cartesian_Grappa_RealTimeCine_Python.xml"
    )
endif()

