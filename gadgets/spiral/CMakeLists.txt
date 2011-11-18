
include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
				       ${TINYXML_INCLUDE_DIR}
				       ${Boost_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetronspiral SHARED SpiralGadget.cpp vds.cpp)

target_link_libraries(gadgetronspiral optimized gpunfft ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES})
target_link_libraries(gadgetronspiral debug gpunfft ${ACE_DEBUG_LIBRARY} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES})

install (TARGETS gadgetronspiral DESTINATION lib)
install (FILES spiral.xml  DESTINATION config)

