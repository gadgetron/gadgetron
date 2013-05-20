IF (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPURADIAL__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/cgsense
  ${CMAKE_SOURCE_DIR}/gadgets/core
  )

add_library(gadgetron_gpuradial SHARED 
#  cuRadialSensPrepGadget.cpp 
  ${ISMRMRD_XSD_SOURCE})

target_link_libraries(gadgetron_gpuradial
  gpunfft gpusolvers gpuoperators gpuparallelmri
  ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpuradial DESTINATION lib)
