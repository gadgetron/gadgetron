IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_RADIAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/pmri
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
  )

include_directories(${CUDA_INCLUDE_DIRS})

if (ARMADILLO_FOUND)
  list(APPEND PHASE_GADGET RadialPhaseCorrectionGadget.cpp)
elseif (ARMADILLO_FOUND)
  MESSAGE("Armadillo not found, not compiling radial phase correction gadget")
endif (ARMADILLO_FOUND)

add_library(gadgetron_gpuradial SHARED 
  gpuRadialPrepGadget.h gpuRadialPrepGadget.cpp 
  gpuRadialSensePrepGadget.h gpuRadialSensePrepGadget.cpp 
  gpuRadialSpiritPrepGadget.h gpuRadialSpiritPrepGadget.cpp 
  gpuRetroGatedSensePrepGadget.h gpuRetroGatedSensePrepGadget.cpp
  ${PHASE_GADGET}
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_gpuradial
  gpunfft gpusolvers gpuoperators gpuparallelmri cpucore cpucore_math gpucore
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpuradial DESTINATION lib)

add_subdirectory(config)
