
include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
				       ${TINYXML_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/solvers
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetroncgsense SHARED     GPUCGGadget.cpp
								GPUCGGoldenRadial.cpp)

target_link_libraries(gadgetroncgsense solvers gpuparallelmri ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES})

install (TARGETS gadgetroncgsense DESTINATION lib)
install (FILES radial_single.xml  DESTINATION config)

