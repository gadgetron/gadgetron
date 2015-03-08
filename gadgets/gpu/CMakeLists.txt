IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUGADGET__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

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

add_library(gadgetron_gpugadget SHARED 
 cuFFTGadget.h cuFFTGadget.cpp) 

set_target_properties(gadgetron_gpugadget PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_gpugadget
  gadgetron_gadgetbase
  gadgetron_toolbox_gpunfft gadgetron_toolbox_gpusolvers gadgetron_toolbox_gpuoperators gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )


install (FILES  cuFFTGadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_gpugadget DESTINATION lib COMPONENT main)

