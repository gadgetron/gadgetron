include_directories(${ACE_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
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
								NoiseAdjustGadget.cpp)

target_link_libraries(gadgetroncore ${ACE_LIBRARIES} ${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES})

install(TARGETS gadgetroncore DESTINATION lib)
install(FILES default.xml DESTINATION config)
