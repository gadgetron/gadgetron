# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/config" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_realtime.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_realtime.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_realtime.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_gpusense_cg.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_gpusense_cg.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_cg.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode3_gpusense_cg.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_gpusense_sb.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_gpusense_sb.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_sb.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_nlcg.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode3_gpusense_sb.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_gpu_ktsense.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_gpu_ktsense.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpu_ktsense.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode3_os_realtime.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/spirit.xml"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/config" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_gpusense_cg_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_gpusense_cg_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_cg_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_nlcg_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode0_gpusense_sb_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/fixed_radial_mode1_gpusense_sb_unoptimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/radial/config/golden_radial_mode2_gpusense_sb_unoptimized.xml"
    )
endif()

