if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUSENSE__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/gadgets/sense
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ISMRMRD_XSD_INCLUDE_DIR}
)

add_library(gadgetron_gpusense SHARED 
  gpuCgSenseGadget.cpp 
  gpuCgKtSenseGadget.cpp 
  gpuSbSenseGadget.cpp 
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_gpusense 
  cpucore gpucore gpusolvers gpuoperators gpuparallelmri 
  ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpusense DESTINATION lib)

install (FILES 
  SenseJob.h
  DESTINATION include)
