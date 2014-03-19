if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUPMRI__)
  # The two flags below is to fix Windows problems in relation to multiple defined constructors in our headers
  SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE") 
  SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE") 
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
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ISMRMRD_XSD_INCLUDE_DIR}
)

add_library(gadgetron_gpuparallelmri SHARED 
  gpuCgSenseGadget.cpp 
  gpuCgKtSenseGadget.cpp 
  gpuSbSenseGadget.cpp 
  gpuGenericSensePrepGadget.cpp
  gpuCgSpiritGadget.cpp 
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_gpuparallelmri 
  cpucore gpucore gpusolvers gpuoperators gpuparallelmri 
  ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_gpuparallelmri DESTINATION lib)

install (FILES 
  DESTINATION include)

add_subdirectory(config)
