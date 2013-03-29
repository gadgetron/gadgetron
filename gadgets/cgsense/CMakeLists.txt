if (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CGSENSE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
)

add_library(gadgetroncgsense SHARED GPUCGGadget.cpp GPUCGGoldenRadial.cpp GPUCGFixedRadial.cpp GPUCGGadgetGeneric.cpp)

target_link_libraries(gadgetroncgsense gpusolvers gpuoperators gpuparallelmri tinyxml ${ISMRMRD_LIBRARIES} ${ISMRMRD_XSD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (TARGETS gadgetroncgsense DESTINATION lib)
install (FILES radial_single.xml  DESTINATION config)
