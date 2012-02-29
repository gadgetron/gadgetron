if (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GRAPPA__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

include_directories(   ${ACE_INCLUDE_DIR}
				       ${CUDA_INCLUDE_DIRS}
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetrongrappa SHARED GrappaGadget.cpp
								   GrappaCalibrationBuffer.cpp
							 	   GrappaWeights.cpp
							  	   GrappaWeightsCalculator.cpp
							  	   GrappaUnmixingGadget.cpp)

target_link_libraries(gadgetrongrappa gpuparallelmri cpucore tinyxml ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

install (TARGETS gadgetrongrappa DESTINATION lib)
install (FILES grappa.xml DESTINATION config)


