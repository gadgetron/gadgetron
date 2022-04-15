if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRICORE__)
endif ()

find_package(ZFP)

include_directories(${HDF5_INCLUDE_DIRS})

if(MSVC)
  set_source_files_properties(CompressedFloatBufferAvx2.cpp PROPERTIES COMPILE_FLAGS "/arch:AVX2")
else()
  set_source_files_properties(CompressedFloatBufferSse41.cpp PROPERTIES COMPILE_FLAGS "-msse4.1")
  set_source_files_properties(CompressedFloatBufferAvx2.cpp PROPERTIES COMPILE_FLAGS "-mavx2 -mfma")
endif()


set(gadgetron_mricore_header_files GadgetMRIHeaders.h
        NoiseAdjustGadget.h
        PCACoilGadget.h
        AccumulatorGadget.h
        FFTGadget.h
        CombineGadget.h
        CropAndCombineGadget.h
        writers/MRIImageWriter.h
        readers/MRIImageReader.h
        NoiseAdjustGadget_unoptimized.h
        ExtractGadget.h
        FloatToFixPointGadget.h
        RemoveROOversamplingGadget.h
        CoilReductionGadget.h
        ScaleGadget.h
        FlowPhaseSubtractionGadget.h
        readers/GadgetIsmrmrdReader.h
        PhysioInterpolationGadget.h
        AsymmetricEchoAdjustROGadget.h
        MaxwellCorrectionGadget.h
        dependencyquery/DependencyQueryGadget.h
        dependencyquery/DependencyQueryWriter.h
        ComplexToFloatGadget.h
        AcquisitionAccumulateTriggerGadget.h
        BucketToBufferGadget.h
        ImageArraySplitGadget.h
        SimpleReconGadget.h
        ImageSortGadget.h
        generic_recon_gadgets/GenericReconBase.h
        generic_recon_gadgets/GenericReconGadget.h
        generic_recon_gadgets/GenericReconCartesianFFTGadget.h
        generic_recon_gadgets/GenericReconCartesianGrappaGadget.h
        generic_recon_gadgets/GenericReconCartesianSpiritGadget.h
        generic_recon_gadgets/GenericReconCartesianNonLinearSpirit2DTGadget.h
        generic_recon_gadgets/GenericReconCartesianReferencePrepGadget.h
        generic_recon_gadgets/GenericReconPartialFourierHandlingGadget.h
        generic_recon_gadgets/GenericReconPartialFourierHandlingFilterGadget.h
        generic_recon_gadgets/GenericReconPartialFourierHandlingPOCSGadget.h
        generic_recon_gadgets/GenericReconKSpaceFilteringGadget.h
        generic_recon_gadgets/GenericReconFieldOfViewAdjustmentGadget.h
        generic_recon_gadgets/GenericReconImageArrayScalingGadget.h
        generic_recon_gadgets/GenericReconEigenChannelGadget.h
        generic_recon_gadgets/GenericReconNoiseStdMapComputingGadget.h
        generic_recon_gadgets/GenericImageReconGadget.h
        generic_recon_gadgets/GenericReconAccumulateImageTriggerGadget.h
        generic_recon_gadgets/GenericImageReconArrayToImageGadget.h
        generic_recon_gadgets/GenericReconReferenceKSpaceDelayedBufferGadget.h
        generic_recon_gadgets/GenericReconImageToImageArrayGadget.h 
        ImageFinishGadget.h
        dependencyquery/NoiseSummaryGadget.h
        NHLBICompression.h
    	  cpuisa.h
        ImageAccumulatorGadget.h
        writers/GadgetIsmrmrdWriter.h
        ImageResizingGadget.h
        ImageArraySendMixin.h
        ImageArraySendMixin.hpp
        DenoiseGadget.h
        CoilComputationGadget.h 
        ImageFFTGadget.h
        FlagTriggerGadget.h
        ImageIndexGadget.h
        IsmrmrdDumpGadget.h
        )

set(gadgetron_mricore_src_files
        NoiseAdjustGadget.cpp
        PCACoilGadget.cpp
        AccumulatorGadget.cpp
        FFTGadget.cpp
        CombineGadget.cpp
        CropAndCombineGadget.cpp
        ImageFinishGadget.cpp
        writers/MRIImageWriter.cpp
        readers/MRIImageReader.cpp
        NoiseAdjustGadget_unoptimized.cpp
        ExtractGadget.cpp
        FloatToFixPointGadget.cpp
        RemoveROOversamplingGadget.cpp
        CoilReductionGadget.cpp
        ScaleGadget.cpp
        FlowPhaseSubtractionGadget.cpp
        readers/GadgetIsmrmrdReader.cpp
        PhysioInterpolationGadget.cpp
        AsymmetricEchoAdjustROGadget.cpp
        MaxwellCorrectionGadget.cpp
        dependencyquery/DependencyQueryGadget.cpp
        dependencyquery/DependencyQueryWriter.cpp
        ComplexToFloatGadget.cpp
        AcquisitionAccumulateTriggerGadget.cpp
        BucketToBufferGadget.cpp
        ImageArraySplitGadget.cpp
        SimpleReconGadget.cpp
        ImageSortGadget.cpp
        generic_recon_gadgets/GenericReconBase.cpp
        generic_recon_gadgets/GenericReconGadget.cpp
        generic_recon_gadgets/GenericReconCartesianFFTGadget.cpp
        generic_recon_gadgets/GenericReconCartesianGrappaGadget.cpp
        generic_recon_gadgets/GenericReconCartesianSpiritGadget.cpp
        generic_recon_gadgets/GenericReconCartesianNonLinearSpirit2DTGadget.cpp
        generic_recon_gadgets/GenericReconCartesianReferencePrepGadget.cpp
        generic_recon_gadgets/GenericReconPartialFourierHandlingGadget.cpp
        generic_recon_gadgets/GenericReconPartialFourierHandlingFilterGadget.cpp
        generic_recon_gadgets/GenericReconPartialFourierHandlingPOCSGadget.cpp
        generic_recon_gadgets/GenericReconReferenceKSpaceDelayedBufferGadget.cpp
        generic_recon_gadgets/GenericReconKSpaceFilteringGadget.cpp
        generic_recon_gadgets/GenericReconFieldOfViewAdjustmentGadget.cpp
        generic_recon_gadgets/GenericReconImageArrayScalingGadget.cpp
        generic_recon_gadgets/GenericReconEigenChannelGadget.cpp
        generic_recon_gadgets/GenericReconNoiseStdMapComputingGadget.cpp
        generic_recon_gadgets/GenericImageReconGadget.cpp
        generic_recon_gadgets/GenericReconAccumulateImageTriggerGadget.cpp
        generic_recon_gadgets/GenericImageReconArrayToImageGadget.cpp
        generic_recon_gadgets/GenericReconImageToImageArrayGadget.cpp 
        CompressedFloatBuffer.cpp
        CompressedFloatBufferSse41.cpp
        CompressedFloatBufferAvx2.cpp
        cpuisa.cpp
        dependencyquery/NoiseSummaryGadget.cpp
        ImageAccumulatorGadget.cpp
        writers/GadgetIsmrmrdWriter.cpp
        ImageResizingGadget.cpp
        AutoScaleGadget.cpp
        DenoiseGadget.cpp
        CoilComputationGadget.cpp
        ImageFFTGadget.cpp
        FlagTriggerGadget.cpp
        ImageIndexGadget.cpp
        IsmrmrdDumpGadget.cpp
        )

set(gadgetron_mricore_config_files
        config/default.xml
        config/default_short.xml
        config/default_optimized.xml
        config/default_measurement_dependencies.xml
        config/default_measurement_dependencies_ismrmrd_storage.xml
        config/isalive.xml
        config/gtquery.xml
        generic_recon_gadgets/config/Generic_Cartesian_FFT.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_SNR.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_T2W.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_RealTimeCine.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_RealTimeCine_Cloud.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_EPI.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_EPI_AVE.xml
        generic_recon_gadgets/config/Generic_Cartesian_Spirit.xml
        generic_recon_gadgets/config/Generic_Cartesian_Spirit_RealTimeCine.xml
        generic_recon_gadgets/config/Generic_Cartesian_Spirit_SASHA.xml
        generic_recon_gadgets/config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine.xml
        generic_recon_gadgets/config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine.xml
        generic_recon_gadgets/config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml
        generic_recon_gadgets/config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine_Cloud.xml
        generic_recon_gadgets/config/Generic_Cartesian_Grappa_Cine_Denoise.xml
        generic_recon_gadgets/config/Generic_Cartesian_Image_Chain_FFT.xml
        config/NoiseSummary.xml
        )

# if there is python, add gadgets using python
if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)

    message("Build python generic gadgets ... ")
    include_directories(
      ${CMAKE_SOURCE_DIR}/core
      ${CMAKE_SOURCE_DIR}/toolboxes/python
      ${PYTHON_INCLUDE_PATH}
      ${NUMPY_INCLUDE_DIRS}
      )

      set( gadgetron_mricore_header_files ${gadgetron_mricore_header_files} generic_recon_gadgets/GenericReconCartesianGrappaAIGadget.h)
      set( gadgetron_mricore_src_files ${gadgetron_mricore_src_files} generic_recon_gadgets/GenericReconCartesianGrappaAIGadget.cpp)
      set( gadgetron_mricore_config_files ${gadgetron_mricore_config_files} generic_recon_gadgets/config/Generic_Cartesian_Grappa_AI.xml)

      set( gadgetron_python_models_files models/grappa_ai.py )
endif ()

add_library(gadgetron_mricore SHARED
    gadgetron_mricore_export.h
    ${gadgetron_mricore_header_files}
    ${gadgetron_mricore_src_files}
    ${gadgetron_mricore_config_files}
    DenoiseGadget.cpp DenoiseGadget.h)

set_target_properties(gadgetron_mricore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mricore
        gadgetron_core
        gadgetron_core_writers
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpufft
        gadgetron_toolbox_image_analyze_io
        gadgetron_toolbox_denoise
        gadgetron_toolbox_hostutils
        gadgetron_toolbox_cpuoperator
        gadgetron_toolbox_cpuklt
        gadgetron_toolbox_mri_core
        )

target_include_directories(gadgetron_mricore
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

if (ZFP_FOUND)
    message("ZFP Found")
    add_definitions(-DGADGETRON_COMPRESSION_ZFP)
    target_include_directories(
            gadgetron_mricore
            PUBLIC
            ${ZFP_INCLUDE_DIR})

    target_link_libraries(gadgetron_mricore ${ZFP_LIBRARIES})
else ()
    message("ZFP NOT Found")
endif ()

if (Boost_PYTHON3_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
   target_link_libraries(gadgetron_mricore
                        gadgetron_toolbox_python
                         )

    install(FILES
        ${gadgetron_python_models_files}
        DESTINATION ${GADGETRON_INSTALL_PYTHON_MODULE_PATH} COMPONENT main)
endif ()

install(FILES
        gadgetron_mricore_export.h
        ${gadgetron_mricore_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${gadgetron_mricore_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_mricore
        EXPORT gadgetron-export
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)
