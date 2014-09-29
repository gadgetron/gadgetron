IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_SPIRAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CUDA_INCLUDE_DIRS}
  )

add_library(gadgetron_spiral SHARED 
  gadgetron_spiral_export.h 
  vds.cpp 
  gpuSpiralSensePrepGadget.h gpuSpiralSensePrepGadget.cpp 
  SpiralToGenericGadget.h SpiralToGenericGadget.cpp)

set_target_properties(gadgetron_spiral PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_spiral
  gadgetron_toolbox_cpucore gadgetron_toolbox_gpucore gadgetron_toolbox_gpunfft gadgetron_toolbox_gpusolvers gadgetron_toolbox_gpuoperators
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

target_link_libraries(gadgetron_spiral gadgetron_toolbox_gpuparallelmri)

install (TARGETS gadgetron_spiral DESTINATION lib)
install (FILES vds.h gadgetron_spiral_export.h 
                     gpuSpiralSensePrepGadget.h 
                     SpiralToGenericGadget.h 
                     DESTINATION include)

add_subdirectory(config)
