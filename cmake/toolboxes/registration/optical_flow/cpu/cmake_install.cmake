# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu

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
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoCKOpticalFlowSolver.cpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoCKOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoHSOpticalFlowSolver.cpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoHSOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoOpticalFlowSolver.cpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoLinearResampleOperator.cpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/hoLinearResampleOperator.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu/libgadgetron_toolbox_cpureg.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu/libgadgetron_toolbox_cpureg.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu/libgadgetron_toolbox_cpureg.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_cpureg.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegParametricTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegHomogenousTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegRigid2DTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegRigid3DTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegNonParametricTransformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/transformation/hoImageRegDeformationField.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegParametricSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegParametricDownHillSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegParametricGradientDescentSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegNonParametricSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegDeformationFieldSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/solver/hoImageRegDeformationFieldBidirectionalSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/warper/hoImageRegWarper.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilarity.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilarityHistogramBased.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilarityLocalCCR.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilarityMutualInformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilarityNormalizedMutualInformation.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/dissimilarity/hoImageRegDissimilaritySSD.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/register/hoImageRegRegister.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/register/hoImageRegParametricRegister.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/register/hoImageRegNonParametricRegister.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/register/hoImageRegDeformationFieldRegister.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/register/hoImageRegDeformationFieldBidirectionalRegister.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/cpu/application/hoImageRegContainer2DRegistration.h"
    )
endif()

