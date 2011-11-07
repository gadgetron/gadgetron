if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CGSENSE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

include_directories(		       ${ACE_INCLUDE_DIR}
				       ${CUDA_INCLUDE_DIRS}
				       ${TINYXML_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
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

target_link_libraries(gadgetroncgsense optimized solvers gpuparallelmri ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES})
target_link_libraries(gadgetroncgsense debug solvers gpuparallelmri ${ACE_DEBUG_LIBRARY} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES})

install (TARGETS gadgetroncgsense DESTINATION lib)
install (FILES radial_single.xml  DESTINATION config)

