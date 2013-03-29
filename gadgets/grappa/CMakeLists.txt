if (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GRAPPA__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
)

add_library(gadgetrongrappa SHARED GrappaGadget.cpp
								   GrappaCalibrationBuffer.cpp
							 	   GrappaWeights.cpp
							  	   GrappaWeightsCalculator.cpp
							  	   GrappaUnmixingGadget.cpp)

target_link_libraries(gadgetrongrappa gpuparallelmri cpucore tinyxml ${ISMRMRD_LIBRARIES}  ${ISMRMRD_XSD_LIBRARIES} ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

install (TARGETS gadgetrongrappa DESTINATION lib)
install (FILES grappa.xml DESTINATION config)


