if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GRAPPA__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
)

add_library(gadgetron_grappa SHARED 
	GrappaGadget.cpp
	GrappaCalibrationBuffer.cpp
	GrappaWeights.cpp
	GrappaWeightsCalculator.cpp
	GrappaUnmixingGadget.cpp
	${ISMRMRD_XSD_SOURCE}
    )

target_link_libraries(gadgetron_grappa 
  gpuparallelmri cpucore
  ${Boost_LIBRARIES}
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${XERCESC_LIBRARIES} 
  )

install (TARGETS gadgetron_grappa DESTINATION lib)
install (FILES grappa.xml grappa_float.xml DESTINATION config)
