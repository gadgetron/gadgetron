IF (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

if (MKL_FOUND)
    # This is a fix for the bug in SVD when MKL is multi-threaded
    MESSAGE("MKL Found, enabling MKL for mri_core gadgets.")
    add_definitions(-DHAVE_MKL)
    # These are needed to get the linking to work properly when
    # MKL is installed, but Armadillo is NOT using it.
    list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
    list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
endif (MKL_FOUND)

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
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
    ${ARMADILLO_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR}
    ${ISMRMRD_INCLUDE_DIR}
)

if (ARMADILLO_FOUND)
    list(APPEND OPTIMIZED_GADGETS NoiseAdjustGadget.cpp)
    list(APPEND OPTIMIZED_GADGETS PCACoilGadget.cpp)
    list(APPEND OPTIMIZED_GADGET_HEADERS NoiseAdjustGadget.h)
    list(APPEND OPTIMIZED_GADGET_HEADERS PCACoilGadget.h)
    list(APPEND OPTIMIZED_GADGET_LIBS gadgetron_toolbox_cpucore_math ${ARMADILLO_LIBRARIES})
elseif (ARMADILLO_FOUND)
    MESSAGE("Armadillo not found, i.e. not compiling Blas/Lapack optimized MRI Gadgets")
endif (ARMADILLO_FOUND)

set( gadgetron_mricore_header_files GadgetMRIHeaders.h 
                                    AcquisitionPassthroughGadget.h 
                                    AcquisitionFinishGadget.h 
                                    AccumulatorGadget.h 
                                    FFTGadget.h 
                                    ImageFinishGadget.h 
                                    CombineGadget.h
                                    CropAndCombineGadget.h  
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
                                    SimpleReconGadget.h
                                    ImageSortGadget.h)

set( gadgetron_mricore_src_files AcquisitionPassthroughGadget.cpp 
                                AcquisitionFinishGadget.cpp 
                                AccumulatorGadget.cpp 
                                FFTGadget.cpp 
                                ImageFinishGadget.cpp 
                                CombineGadget.cpp 
                                CropAndCombineGadget.cpp 
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
                                SimpleReconGadget.cpp
                                ImageSortGadget.cpp)

if (WIN32)
    set( gadgetron_mricore_header_files ${gadgetron_mricore_header_files} WhiteNoiseInjectorGadget.h )
    set( gadgetron_mricore_src_files ${gadgetron_mricore_src_files} WhiteNoiseInjectorGadget.cpp )
endif (WIN32)

set( gadgetron_mricore_config_files
    default.xml
    default_short.xml
    default_optimized.xml
    default_measurement_dependencies.xml
)

add_library(gadgetron_mricore SHARED 
    gadgetron_mricore_export.h 
    ${gadgetron_mricore_header_files} 
    ${gadgetron_mricore_src_files}
    ${gadgetron_mricore_config_files}
    ${OPTIMIZED_GADGETS}
    ${OPTIMIZED_GADGET_HEADERS}
)

set_target_properties(gadgetron_mricore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_mricore 
    gadgetron_gadgetbase
    gadgetron_toolbox_log
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuklt
    ${ISMRMRD_LIBRARIES} 
    ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${BOOST_LIBRARIES}
    ${OPTIMIZED_GADGET_LIBS}
)

install(FILES 
    gadgetron_mricore_export.h
    ${gadgetron_mricore_header_files}
    ${OPTIMIZED_GADGET_HEADERS}
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

if (ARMADILLO_FOUND)
    install(FILES ${gadgetron_mricore_config_files} DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)
endif (ARMADILLO_FOUND)

install(FILES ismrmrd_dump.xml DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT main)

install(TARGETS gadgetron_mricore DESTINATION lib COMPONENT main)
