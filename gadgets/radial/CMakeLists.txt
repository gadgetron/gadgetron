IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_RADIAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/gadgets/sense
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${ISMRMRD_XSD_INCLUDE_DIR}
  )

if(CUDA_FOUND)
  include_directories(${CUDA_INCLUDE_DIRS})
  
  add_library(gadgetron_gpuradial SHARED 
    gpuRadialSensePrepGadget.cpp 
    gpuRetroGatedSensePrepGadget.cpp
    ${ISMRMRD_XSD_SOURCE})
  
  target_link_libraries(gadgetron_gpuradial
    gpunfft gpusolvers gpuoperators gpuparallelmri cpucore cpucore_math gpucore
    ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
    optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
    )
  
  install (TARGETS gadgetron_gpuradial DESTINATION lib)
endif(CUDA_FOUND)

add_subdirectory(config)
