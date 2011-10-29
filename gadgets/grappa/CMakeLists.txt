
include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
				       ${TINYXML_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetrongrappa SHARED     GrappaGadget.cpp
								GrappaCalibrationBuffer.cpp
							 	GrappaWeights.cpp
							  	GrappaWeightsCalculator.cpp
								RemoveROOversamplingGadget.cpp)

target_link_libraries(gadgetrongrappa gpuparallelmri ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES})

install (TARGETS gadgetrongrappa DESTINATION lib)
install (FILES grappa.xml DESTINATION config)


