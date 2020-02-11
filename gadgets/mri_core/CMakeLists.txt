if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRICORE__)
endif ()

find_package(ZFP)




set(gadgetron_mricore_header_files GadgetMRIHeaders.h
        NoiseAdjustGadget.h
        PCACoilGadget.h
        RateLimitGadget.h
        AcquisitionPassthroughGadget.h
        AccumulatorGadget.h
        FFTGadget.h
        CombineGadget.h
        CropAndCombineGadget.h
        ImageWriterGadget.h
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
        IsmrmrdDumpGadget.h
        AsymmetricEchoAdjustROGadget.h
        MaxwellCorrectionGadget.h
        CplxDumpGadget.h
        dependencyquery/DependencyQueryGadget.h
        dependencyquery/DependencyQueryWriter.h
        ComplexToFloatGadget.h
        AcquisitionAccumulateTriggerGadget.h
        BucketToBufferGadget.h
        ImageArraySplitGadget.h
        PseudoReplicatorGadget.h
        SimpleReconGadget.h
        ImageSortGadget.h
        GenericReconBase.h
        GenericReconGadget.h
        GenericReconCartesianFFTGadget.h
        GenericReconCartesianGrappaGadget.h
        GenericReconCartesianSpiritGadget.h
        GenericReconCartesianNonLinearSpirit2DTGadget.h
        GenericReconCartesianReferencePrepGadget.h
        GenericReconPartialFourierHandlingGadget.h
        GenericReconPartialFourierHandlingFilterGadget.h
        GenericReconPartialFourierHandlingPOCSGadget.h
        GenericReconKSpaceFilteringGadget.h
        GenericReconFieldOfViewAdjustmentGadget.h
        GenericReconImageArrayScalingGadget.h
        GenericReconEigenChannelGadget.h
        GenericReconNoiseStdMapComputingGadget.h
        GenericImageReconGadget.h
        GenericReconAccumulateImageTriggerGadget.h
        GenericImageReconArrayToImageGadget.h
        GenericReconReferenceKSpaceDelayedBufferGadget.h
        WhiteNoiseInjectorGadget.h
        ImageFinishGadget.h
        dependencyquery/NoiseSummaryGadget.h
        NHLBICompression.h
        ImageAccumulatorGadget.h
        writers/GadgetIsmrmrdWriter.h
        ImageResizingGadget.h
        DenoiseGadget.h)

set(gadgetron_mricore_src_files
        AcquisitionPassthroughGadget.cpp
        NoiseAdjustGadget.cpp
        PCACoilGadget.cpp
        AccumulatorGadget.cpp
        FFTGadget.cpp
        CombineGadget.cpp
        CropAndCombineGadget.cpp
        ImageFinishGadget.cpp
        ImageWriterGadget.cpp
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
        IsmrmrdDumpGadget.cpp
        AsymmetricEchoAdjustROGadget.cpp
        MaxwellCorrectionGadget.cpp
        CplxDumpGadget.cpp
        dependencyquery/DependencyQueryGadget.cpp
        dependencyquery/DependencyQueryWriter.cpp
        ComplexToFloatGadget.cpp
        AcquisitionAccumulateTriggerGadget.cpp
        BucketToBufferGadget.cpp
        ImageArraySplitGadget.cpp
        PseudoReplicatorGadget.cpp
        SimpleReconGadget.cpp
        ImageSortGadget.cpp
        GenericReconBase.cpp
        GenericReconGadget.cpp
        GenericReconCartesianFFTGadget.cpp
        GenericReconCartesianGrappaGadget.cpp
        GenericReconCartesianSpiritGadget.cpp
        GenericReconCartesianNonLinearSpirit2DTGadget.cpp
        GenericReconCartesianReferencePrepGadget.cpp
        GenericReconPartialFourierHandlingGadget.cpp
        GenericReconPartialFourierHandlingFilterGadget.cpp
        GenericReconPartialFourierHandlingPOCSGadget.cpp
        GenericReconReferenceKSpaceDelayedBufferGadget.cpp
        GenericReconKSpaceFilteringGadget.cpp
        GenericReconFieldOfViewAdjustmentGadget.cpp
        GenericReconImageArrayScalingGadget.cpp
        GenericReconEigenChannelGadget.cpp
        GenericReconNoiseStdMapComputingGadget.cpp
        GenericImageReconGadget.cpp
        GenericReconAccumulateImageTriggerGadget.cpp
        GenericImageReconArrayToImageGadget.cpp
        WhiteNoiseInjectorGadget.cpp
        RateLimitGadget.cpp
        dependencyquery/NoiseSummaryGadget.cpp
        ImageAccumulatorGadget.cpp
        writers/GadgetIsmrmrdWriter.cpp
        ImageResizingGadget.cpp
        AutoScaleGadget.cpp
        DenoiseGadget.cpp)

set(gadgetron_mricore_config_files
        config/default.xml
        config/default_short.xml
        config/default_optimized.xml
        config/default_measurement_dependencies.xml
        config/default_measurement_dependencies_ismrmrd_storage.xml
        config/isalive.xml
        config/gtquery.xml
        config/Generic_Cartesian_FFT.xml
        config/Generic_Cartesian_Grappa.xml
        config/Generic_Cartesian_Grappa_SNR.xml
        config/Generic_Cartesian_Grappa_T2W.xml
        config/Generic_Cartesian_Grappa_RealTimeCine.xml
        config/Generic_Cartesian_Grappa_RealTimeCine_Cloud.xml
        config/Generic_Cartesian_Grappa_EPI.xml
        config/Generic_Cartesian_Grappa_EPI_AVE.xml
        config/Generic_Cartesian_Spirit.xml
        config/Generic_Cartesian_Spirit_RealTimeCine.xml
        config/Generic_Cartesian_Spirit_SASHA.xml
        config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine.xml
        config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine.xml
        config/Generic_Cartesian_NonLinear_Spirit_RealTimeCine_Cloud.xml
        config/Generic_Cartesian_RandomSampling_NonLinear_Spirit_RealTimeCine_Cloud.xml
        config/Generic_Cartesian_Grappa_Cine_Denoise.xml
        config/NoiseSummary.xml
        config/ismrmrd_dump.xml
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

      set( gadgetron_mricore_header_files ${gadgetron_mricore_header_files} GenericReconCartesianGrappaAIGadget.h)
      set( gadgetron_mricore_src_files ${gadgetron_mricore_src_files} GenericReconCartesianGrappaAIGadget.cpp)
      set( gadgetron_mricore_config_files ${gadgetron_mricore_config_files} config/Generic_Cartesian_Grappa_AI.xml)

endif ()

add_library(gadgetron_mricore SHARED
    gadgetron_mricore_export.h
    ${gadgetron_mricore_header_files}
    ${gadgetron_mricore_src_files}
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
        boost
        ISMRMRD
        )

target_include_directories(gadgetron_mricore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

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
endif ()


install(FILES
        gadgetron_mricore_export.h
        ${gadgetron_mricore_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(FILES ${gadgetron_mricore_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_mricore
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


