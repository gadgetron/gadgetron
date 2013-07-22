if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_MOCO__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   
  ${CMAKE_SOURCE_DIR}/gadgets/moco
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
#  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
#  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/gpu
#  ${CMAKE_SOURCE_DIR}/toolboxes/operators
#  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ISMRMRD_XSD_INCLUDE_DIR}
)

add_library(gadgetron_moco SHARED
  RegistrationAveragingGadget.cpp
  cpuRegistrationAveragingGadget.cpp
  gpuRegistrationAveragingGadget.cpp
  ${ISMRMRD_XSD_SOURCE}
  )

target_link_libraries(gadgetron_moco
  cpucore gpucore gpusolvers gpuoperators gpuparallelmri 
  ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetron_moco DESTINATION lib)

add_subdirectory(config)
