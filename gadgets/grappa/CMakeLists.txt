if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_GRAPPA__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${ARMADILLO_INCLUDE_DIRS}
)

add_library(gadgetron_grappa SHARED 
    gadgetron_grappa_export.h
    GrappaCalibrationBuffer.h
    GrappaGadget.h
    GrappaUnmixingGadget.h
    GrappaWeights.h
    GrappaWeightsCalculator.h
    GrappaGadget.cpp
    GrappaCalibrationBuffer.cpp
    GrappaWeights.cpp
    GrappaWeightsCalculator.cpp
    GrappaUnmixingGadget.cpp
    ${ISMRMRD_XSD_SOURCE}
    )

target_link_libraries(gadgetron_grappa 
    gpuparallelmri 
    cpucore
    cpufft
    gpufft
    ${Boost_LIBRARIES}
    ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
    ${XERCESC_LIBRARIES} 
  )

install (FILES  gadgetron_grappa_export.h
                GrappaCalibrationBuffer.h
                GrappaGadget.h
                GrappaUnmixingGadget.h
                GrappaWeights.h
                GrappaWeightsCalculator.h 
                DESTINATION include)

install (TARGETS gadgetron_grappa DESTINATION lib)

add_subdirectory(config)
