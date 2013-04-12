if (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GRAPPA__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(Ismrmrd REQUIRED)
find_package(XSD REQUIRED)
find_package(XercesC REQUIRED)

#Process the XSD files
SET(XSDS ${ISMRMRD_SCHEMA_DIR}/ismrmrd.xsd)
SET(XSD_ARGS cxx-tree --generate-serialization)
WRAP_XSD(XSDS_SOURCES XSD_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/schema ${XSDS} OPTIONS ${XSD_ARGS})
INCLUDE_DIRECTORIES(${XSD_INCLUDE_DIR} ${XERCESC_INCLUDE_DIR})

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/gadgets/core
  INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})
)

add_library(gadgetrongrappa SHARED 
	GrappaGadget.cpp
	GrappaCalibrationBuffer.cpp
	GrappaWeights.cpp
	GrappaWeightsCalculator.cpp
	GrappaUnmixingGadget.cpp
	${XSDS_SOURCES}
    )

target_link_libraries(gadgetrongrappa 
  gpuparallelmri cpucore tinyxml 
  ${Boost_LIBRARIES}
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  ${XERCESC_LIBRARIES} 
  )

install (TARGETS gadgetrongrappa DESTINATION lib)
install (FILES grappa.xml DESTINATION config)


