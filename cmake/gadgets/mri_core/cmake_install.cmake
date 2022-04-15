# Install script for directory: /home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core

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
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/gadgetron_mricore_export.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GadgetMRIHeaders.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/NoiseAdjustGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/PCACoilGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/RateLimitGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/AcquisitionPassthroughGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/AccumulatorGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/FFTGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/CombineGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/CropAndCombineGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageWriterGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/writers/MRIImageWriter.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/readers/MRIImageReader.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/NoiseAdjustGadget_unoptimized.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ExtractGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/FloatToFixPointGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/RemoveROOversamplingGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/CoilReductionGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ScaleGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/FlowPhaseSubtractionGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/readers/GadgetIsmrmrdReader.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/PhysioInterpolationGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/IsmrmrdDumpGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/AsymmetricEchoAdjustROGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/MaxwellCorrectionGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/CplxDumpGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/dependencyquery/DependencyQueryGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/dependencyquery/DependencyQueryWriter.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ComplexToFloatGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/AcquisitionAccumulateTriggerGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/BucketToBufferGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageArraySplitGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/PseudoReplicatorGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/SimpleReconGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageSortGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconBase.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconCartesianFFTGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconCartesianGrappaGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconCartesianSpiritGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconCartesianNonLinearSpirit2DTGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconCartesianReferencePrepGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconPartialFourierHandlingGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconPartialFourierHandlingFilterGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconPartialFourierHandlingPOCSGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconKSpaceFilteringGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconFieldOfViewAdjustmentGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconImageArrayScalingGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconEigenChannelGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconNoiseStdMapComputingGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericImageReconGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconAccumulateImageTriggerGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericImageReconArrayToImageGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconReferenceKSpaceDelayedBufferGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/GenericReconImageToImageArrayGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/WhiteNoiseInjectorGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageFinishGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/dependencyquery/NoiseSummaryGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/NHLBICompression.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/cpuisa.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageAccumulatorGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/writers/GadgetIsmrmrdWriter.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageResizingGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageArraySendMixin.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageArraySendMixin.hpp"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/DenoiseGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/CoilComputationGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageFFTGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/FlagTriggerGadget.h"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/ImageIndexGadget.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmainx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/gadgetron/config" TYPE FILE FILES
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/default.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/default_short.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/default_optimized.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/default_measurement_dependencies.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/default_measurement_dependencies_ismrmrd_storage.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/isalive.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/gtquery.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_FFT.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_SNR.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_T2W.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_RealTimeCine.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_RealTimeCine_Cloud.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_EPI.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_EPI_AVE.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Spirit.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Spirit_RealTimeCine.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Spirit_SASHA.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine_Cloud.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Grappa_Cine_Denoise.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/Generic_Cartesian_Image_Chain_FFT.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/NoiseSummary.xml"
    "/home/javeda2/mrprogs/gadgetron_fork/gadgets/mri_core/config/ismrmrd_dump.xml"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core/libgadgetron_mricore.so.4.2.1"
    "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core/libgadgetron_mricore.so.4.2"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so.4.2.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so.4.2"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/usr/local/lib:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
           NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so"
         RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/javeda2/mrprogs/gadgetron_fork/cmake/gadgets/mri_core/libgadgetron_mricore.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so"
         OLD_RPATH "/home/javeda2/mrprogs/gadgetron_fork/GTBLAS:/home/javeda2/mrprogs/gadgetron_fork/BEFORE:/home/javeda2/mrprogs/gadgetron_fork/INTERFACE:/lib/intel64:/home/javeda2/mrprogs/gadgetron_fork/cmake/core/writers:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/denoise:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/hostutils:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/operators/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/mri_core:/usr/local/lib:/home/javeda2/mrprogs/gadgetron_fork/cmake/core:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/fft/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/klt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/registration/optical_flow/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/image_io:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/dwt/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu/math:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/core/cpu:/home/javeda2/mrprogs/gadgetron_fork/cmake/toolboxes/log:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu:"
         NEW_RPATH ".:/home/javeda2/anaconda3/envs/gadgetron/lib:/lib/intel64:/usr/local/lib:/home/javeda2/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:/usr/lib/x86_64-linux-gnu")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/javeda2/anaconda3/envs/gadgetron/bin/x86_64-conda-linux-gnu-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libgadgetron_mricore.so")
    endif()
  endif()
endif()

