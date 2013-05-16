if (WIN32)
  ADD_DEFINITIONS(-DTIXML_USE_STL)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_SBSENSE__)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
  ${CMAKE_SOURCE_DIR}/gadgets/cgsense
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
)

add_library(gadgetronsbsense SHARED 
cuSbGadgetGeneric.cpp ${ISMRMRD_XSD_SOURCE})
#  cuSbGadget.cpp cuSbGoldenRadial.cpp cuSbFixedRadial.cpp 

target_link_libraries(gadgetronsbsense 
  cpucore gpucore gpusolvers gpuoperators gpuparallelmri 
  ${Boost_LIBRARIES} ${ISMRMRD_LIBRARIES} ${XERCESC_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

install (TARGETS gadgetronsbsense DESTINATION lib)
install (FILES radial_single_sb.xml DESTINATION config)
