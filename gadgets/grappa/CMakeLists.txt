if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_GRAPPA__)
endif ()

include_directories(
    ${CMAKE_SOURCE_DIR}/gadgets/mri_core
    ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
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
    )

set_target_properties(gadgetron_grappa PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

if(CUDA_FOUND)
    target_link_libraries(gadgetron_grappa 
      gadgetron_gadgetbase
      gadgetron_toolbox_log
      gadgetron_toolbox_gpuparallelmri
      gadgetron_toolbox_cpucore
      gadgetron_toolbox_cpucore_math
      gadgetron_toolbox_cpufft
      gadgetron_toolbox_gpufft
      gadgetron_toolbox_mri_core
      ${Boost_LIBRARIES}
      ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
      optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
      )
else ()
    target_link_libraries(gadgetron_grappa 
      gadgetron_gadgetbase
      gadgetron_toolbox_log
      gadgetron_toolbox_cpucore
      gadgetron_toolbox_cpucore_math
      gadgetron_toolbox_cpufft
      ${Boost_LIBRARIES}
      ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
      optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
      )
    if (ARMADILLO_FOUND)
      target_link_libraries(gadgetron_grappa 
        gadgetron_toolbox_mri_core
      )
     endif()
endif ()

install (FILES  gadgetron_grappa_export.h
                GrappaCalibrationBuffer.h
                GrappaGadget.h
                GrappaUnmixingGadget.h
                GrappaWeights.h
                GrappaWeightsCalculator.h 
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_grappa DESTINATION lib COMPONENT main)

add_subdirectory(config)
