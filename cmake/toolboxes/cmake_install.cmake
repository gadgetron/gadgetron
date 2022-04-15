# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/toolboxes

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
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/libgadgetron_toolbox_gpu.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/libgadgetron_toolbox_gpu.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/libgadgetron_toolbox_gpu.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_toolbox_gpu.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/gadgetron" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/check_CUDA.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/CUBLASContextProvider.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cudaDeviceManager.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_blas.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_converter.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_elemwise.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_fileio.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_math.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_operators.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_reductions.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray_utils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuNDArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/cuSparseMatrix.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/GadgetronCuException.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/GPUTimer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/hoCuNDArray_blas.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/hoCuNDArray_elemwise.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/hoCuNDArray_math.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/hoCuNDArray_utils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/hoCuNDArray.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/radial_utilities.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/real_utilities_device.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/core/gpu/setup_grid.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/dwt/gpu/cuDWTOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/dwt/gpu/cuNDDWT.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/cuFFTCachedPlan.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/cuFFTCachedPlan.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/cuFFTPlan.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/cuFFTPlan.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/cuNDFFT.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/fft/gpu/gpufft_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/hyper/CSI_utils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/hyper/CSIOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/b1_map.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuBuffer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuCartesianSenseOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuNonCartesianKtSenseOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuNonCartesianSenseOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuSenseBuffer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuSenseBufferCg.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuSenseOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuSpiritBuffer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/cuSpiritOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/gpupmri_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/htgrappa.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/sense_utilities.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/pmri/gpu/senseOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/mri/sdc/gpu/cuSDC.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/nfft/gpu/cuGriddingConvolution.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/nfft/gpu/cuNFFT.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuConvolutionOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuDiagonalOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuDiagonalSumOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuDownsampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuFFTOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuImageOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuLaplaceOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuPartialDerivativeOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuPartialDerivativeOperator2.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuTv1dOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuTvOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuTvPicsOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/cuUpsampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/gpuoperators_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuEncodingOperatorContainer.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuIdentityOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuPartialDerivativeOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuTvOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/operators/gpu/hoCuTvPicsOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuCGHSOFSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuCKOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuHSOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuLinearResampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuOpticalFlowSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/cuResampleOperator.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/registration/optical_flow/gpu/gpureg_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuCgPreconditioner.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuGpBbSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuLbfgsSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuLwSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuNlcgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuSbcCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuSbCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuSbcLwSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuSbLwSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/cuSolverUtils.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/gpusolvers_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/hoCuCgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/hoCuGpBbSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/hoCuNlcgSolver.h"
    "/home/javeda2/mrprogs/gadgetron_fork/toolboxes/solvers/gpu/hoCuSbcCgSolver.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/solvers/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/cmr/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_image/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fatwater/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/deblur/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/nfft/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/ffd/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/pattern_recognition/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/T1/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/python/cmake_install.cmake")
  include("/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/plplot/cmake_install.cmake")

endif()

