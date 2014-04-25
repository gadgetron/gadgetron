IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_SPIRAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

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
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  )

add_library(gadgetron_spiral SHARED 
  gadgetron_spiral_export.h 
  vds.cpp 
  gpuSpiralSensePrepGadget.h gpuSpiralSensePrepGadget.cpp 
  SpiralToGenericGadget.h SpiralToGenericGadget.cpp
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_spiral
  cpucore gpucore gpunfft gpusolvers gpuoperators gpuparallelmri
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_spiral DESTINATION lib)
install (FILES vds.h gadgetron_spiral_export.h 
                     gpuSpiralSensePrepGadget.h 
                     SpiralToGenericGadget.h 
                     DESTINATION include)

add_subdirectory(config)
