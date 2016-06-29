if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_RADIAL__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${ARMADILLO_INCLUDE_DIRS}
  )

include_directories(${CUDA_INCLUDE_DIRS})

if (ARMADILLO_FOUND)
  list(APPEND PHASE_GADGET RadialPhaseCorrectionGadget.h RadialPhaseCorrectionGadget.cpp)
  install (FILES  RadialPhaseCorrectionGadget.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
elseif (ARMADILLO_FOUND)
  message("Armadillo not found, not compiling radial phase correction gadget")
endif ()

add_library(gadgetron_gpuradial SHARED 
  gadgetron_radial_export.h
  gpuRadialPrepGadget.h gpuRadialPrepGadget.cpp 
  gpuRadialSensePrepGadget.h gpuRadialSensePrepGadget.cpp 
  gpuRadialSpiritPrepGadget.h gpuRadialSpiritPrepGadget.cpp 
  gpuRetroGatedSensePrepGadget.h gpuRetroGatedSensePrepGadget.cpp
  ${PHASE_GADGET})

set_target_properties(gadgetron_gpuradial PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_gpuradial
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  gadgetron_toolbox_gpunfft gadgetron_toolbox_gpusolvers gadgetron_toolbox_gpuoperators gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

target_link_libraries(gadgetron_gpuradial gadgetron_toolbox_gpuparallelmri )

install (FILES  gadgetron_radial_export.h
                gpuRadialPrepGadget.h 
                gpuRadialSensePrepGadget.h 
                gpuRadialSpiritPrepGadget.h 
                gpuRetroGatedSensePrepGadget.h 
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_gpuradial DESTINATION lib COMPONENT main)

add_subdirectory(config)
