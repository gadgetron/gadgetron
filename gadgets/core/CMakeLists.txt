IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CORE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_library(gadgetroncore SHARED       	AcquisitionPassthroughGadget.cpp
								AcquisitionFinishGadget.cpp 
								AccumulatorGadget.cpp
								FFTGadget.cpp
								ImageFinishGadget.cpp
								CropAndCombineGadget.cpp
								ImageWriterGadget.cpp
								MRIAcquisitionReader.cpp
								MRIImageWriter.cpp
								NoiseAdjustGadget.cpp
								ExtractMagnitudeGadget.cpp)

target_link_libraries(gadgetroncore optimized ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES})
target_link_libraries(gadgetroncore debug ${ACE_DEBUG_LIBRARY} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES})

install(TARGETS gadgetroncore DESTINATION lib)
install(FILES default.xml DESTINATION config)
