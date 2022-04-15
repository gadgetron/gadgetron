# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu

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
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/cpusolver_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoGdSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoCgPreconditioner.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoLsqrSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoGpBbSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoSbCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/hoSolverUtils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/curveFittingSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/HybridLM.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/cpu/simplexLagariaSolver.h"
    )
endif()

