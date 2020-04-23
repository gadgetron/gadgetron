if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MRICORE__)
endif ()

find_package(ZFP)
if (ZFP_FOUND)
    message("ZFP Found")
    add_definitions(-DGADGETRON_COMPRESSION_ZFP)
    include_directories(${ZFP_INCLUDE_DIR})
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
    include_directories( ${MKL_INCLUDE_DIR} )
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif ()

find_package(HDF5 1.8 REQUIRED COMPONENTS C)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/log
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
    ${CMAKE_SOURCE_DIR}/toolboxes/denoise
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
    ${HDF5_C_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
)   

if (ARMADILLO_FOUND)
    list(APPEND OPTIMIZED_GADGETS NoiseAdjustGadget.cpp)
    list(APPEND OPTIMIZED_GADGETS PCACoilGadget.cpp)
    list(APPEND OPTIMIZED_GADGET_HEADERS NoiseAdjustGadget.h)
    list(APPEND OPTIMIZED_GADGET_HEADERS PCACoilGadget.h)
    list(APPEND OPTIMIZED_GADGET_LIBS gadgetron_toolbox_cpucore_math ${ARMADILLO_LIBRARIES})
elseif (ARMADILLO_FOUND)
    message("Armadillo not found, i.e. not compiling Blas/Lapack optimized MRI Gadgets")
endif ()

set( gadgetron_mricore_header_files GadgetMRIHeaders.h
                                    RateLimitGadget.h
                                    AcquisitionPassthroughGadget.h 
                                    AcquisitionFinishGadget.h 
                                    AccumulatorGadget.h 
                                    FFTGadget.h 
                                    ImageFinishGadget.h 
                                    CombineGadget.h
                                    CropAndCombineGadget.h 
                                    CoilSenMapGadget.h
                                    DependencyDataSenderGadget.h
                                    ImageWriterGadget.h 
                                    MRIImageWriter.h
                                    MRIImageReader.h
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
                                    CplxDumpGadget.h 
                                    DependencyQueryGadget.h 
                                    DependencyQueryWriter.h 
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
                                    GenericReconReferenceKSpaceDelayedBufferGadget.h 
                                    GenericImageReconGadget.h 
                                    GenericReconAccumulateImageTriggerGadget.h 
                                    GenericImageReconArrayToImageGadget.h 
                                    WhiteNoiseInjectorGadget.h
                                    NoiseSummaryGadget.h 
                                    NHLBICompression.h
                                    ImageAccumulatorGadget.h
                                    ImageArraySendMixin.h
                                    ImageArraySendMixin.hpp
                                    CoilComputationGadget.h
				    )

set( gadgetron_mricore_src_files AcquisitionPassthroughGadget.cpp 
                                AcquisitionFinishGadget.cpp 
                                AccumulatorGadget.cpp 
                                FFTGadget.cpp 
                                ImageFinishGadget.cpp 
                                CombineGadget.cpp 
                                CropAndCombineGadget.cpp 
                                CoilSenMapGadget.cpp 
                                DependencyDataSenderGadget.cpp 
                                ImageWriterGadget.cpp 
                                MRIImageWriter.cpp
                                MRIImageReader.cpp
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
                                CplxDumpGadget.cpp 
                                DependencyQueryGadget.cpp 
                                DependencyQueryWriter.cpp 
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
                                GenericReconReferenceKSpaceDelayedBufferGadget.cpp 
                                GenericImageReconArrayToImageGadget.cpp 
                                WhiteNoiseInjectorGadget.cpp
                                RateLimitGadget.cpp
                                NoiseSummaryGadget.cpp
                                ImageAccumulatorGadget.cpp
				CoilComputationGadget.cpp
        )

set( gadgetron_mricore_config_files
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
    ${gadgetron_mricore_config_files}
    ${OPTIMIZED_GADGETS}
    ${OPTIMIZED_GADGET_HEADERS}
        DenoiseGadget.cpp DenoiseGadget.h)
# Above: Whitespace misaligned, and one file on each line please.

set_target_properties(gadgetron_mricore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mricore 
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_image_analyze_io
    gadgetron_toolbox_denoise
    gadgetron_toolbox_hostutils
    ${ISMRMRD_LIBRARIES} 
    ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${OPTIMIZED_GADGET_LIBS}
    ${Boost_LIBRARIES}
)
if(ARMADILLO_FOUND)
   target_link_libraries(gadgetron_mricore
    gadgetron_toolbox_cpuoperator
    gadgetron_toolbox_cpuklt
    gadgetron_toolbox_mri_core
)
endif()

if (ZFP_FOUND)
   target_link_libraries(gadgetron_mricore ${ZFP_LIBRARIES})
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
