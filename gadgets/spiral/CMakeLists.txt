
include_directories(   ${ACE_INCLUDE_DIR}
					   ${CUDA_INCLUDE_DIRS}
				       ${Boost_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetronspiral SHARED SpiralGadget.cpp vds.cpp)

target_link_libraries(gadgetronspiral gpunfft tinyxml ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (TARGETS gadgetronspiral DESTINATION lib)
install (FILES spiral.xml  DESTINATION config)

