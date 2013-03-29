find_package(Ismrmrd REQUIRED)

include_directories(
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

add_library(gadgetronspiral SHARED SpiralGadget.cpp SpiralGadgetSW.cpp FlowPhaseSubtractionGadget.cpp vds.cpp)

target_link_libraries(gadgetronspiral gpunfft tinyxml gpusolvers gpuoperators gpuparallelmri ${ISMRMRD_LIBRARIES} ${ISMRMRD_XSD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (TARGETS gadgetronspiral DESTINATION lib)
install (FILES spiral.xml  spiral_MSH.xml DESTINATION config)

