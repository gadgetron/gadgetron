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
        CoilSenMapGadget.h
#        ImageWriterGadget.h
#        MRIImageWriter.h
#        MRIImageReader.h
        NoiseAdjustGadget_unoptimized.h
        ExtractGadget.h
        FloatToFixPointGadget.h
        RemoveROOversamplingGadget.h
        CoilReductionGadget.h
        AutoScaleGadget.h
       FlowPhaseSubtractionGadget.h
        GadgetIsmrmrdReadWrite.h
        PhysioInterpolationGadget.h
        IsmrmrdDumpGadget.h
        AsymmetricEchoAdjustROGadget.h
        MaxwellCorrectionGadget.h
#        CplxDumpGadget.h
#        DependencyQueryGadget.h
#        DependencyQueryWriter.h
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
        WhiteNoiseInjectorGadget.h
#        NoiseSummaryGadget.h
        NHLBICompression.h
        ImageAccumulatorGadget.h
        DenoiseGadget.h)

set(gadgetron_mricore_src_files
        AcquisitionPassthroughGadget.cpp
        NoiseAdjustGadget.cpp
        PCACoilGadget.cpp
        AccumulatorGadget.cpp
        FFTGadget.cpp
        CombineGadget.cpp
        CropAndCombineGadget.cpp
        CoilSenMapGadget.cpp
#        ImageWriterGadget.cpp
#        MRIImageWriter.cpp
#        MRIImageReader.cpp
        NoiseAdjustGadget_unoptimized.cpp
        ExtractGadget.cpp
        FloatToFixPointGadget.cpp
        RemoveROOversamplingGadget.cpp
        CoilReductionGadget.cpp
        AutoScaleGadget.cpp
        FlowPhaseSubtractionGadget.cpp
        GadgetIsmrmrdReadWrite.cpp
        PhysioInterpolationGadget.cpp
        IsmrmrdDumpGadget.cpp
        AsymmetricEchoAdjustROGadget.cpp
        MaxwellCorrectionGadget.cpp
#        CplxDumpGadget.cpp
#        DependencyQueryGadget.cpp
#        DependencyQueryWriter.cpp
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
#        NoiseSummaryGadget.cpp
        ImageAccumulatorGadget.cpp
        DenoiseGadget.cpp)

set(gadgetron_mricore_config_files
        config/default.xml
        config/default_short.xml
        config/default_optimized.xml
        config/default_measurement_dependencies.xml
        config/default_measurement_dependencies_ismrmrd_storage.xml
        config/default_measurement_dependencies_Noise_CoilSen_SCC.xml
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
        config/NoiseSummary.xml
        config/ismrmrd_dump.xml
        )

add_library(gadgetron_mricore SHARED
        gadgetron_mricore_export.h
        ${gadgetron_mricore_header_files}
        ${gadgetron_mricore_src_files}
        ${gadgetron_mricore_config_files})

set_target_properties(gadgetron_mricore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mricore
        gadgetron_core
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpufft
        gadgetron_toolbox_image_analyze_io
        gadgetron_toolbox_denoise
        gadgetron_toolbox_hostutils
        gadgetron_toolbox_cpuoperator
        gadgetron_toolbox_cpuklt
        gadgetron_toolbox_mri_core
        ${ISMRMRD_LIBRARIES}
        ${FFTW3_LIBRARIES}
        boost
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


if (MKL_FOUND)
    # This is a fix for the bug in SVD when MKL is multi-threaded
    message("MKL Found, enabling MKL for mri_core gadgets.")
    add_definitions(-DHAVE_MKL)
    # These are needed to get the linking to work properly when
    # MKL is installed, but Armadillo is NOT using it.
    list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
    list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
    target_include_directories(gadgetron_mricore
            PUBLIC
            ${MKL_INCLUDE_DIR})
    link_directories(${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR})
endif ()

install(FILES
        gadgetron_mricore_export.h
        ${gadgetron_mricore_header_files}
        ${OPTIMIZED_GADGET_HEADERS}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if (ARMADILLO_FOUND)
    install(FILES ${gadgetron_mricore_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
endif ()

install(TARGETS gadgetron_mricore DESTINATION lib COMPONENT main)
