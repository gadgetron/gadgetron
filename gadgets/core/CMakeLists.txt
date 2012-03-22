IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CORE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
					   ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools)

add_library(gadgetroncore SHARED AcquisitionPassthroughGadget.cpp
								 AcquisitionFinishGadget.cpp 
								 AccumulatorGadget.cpp
								 FFTGadget.cpp
								 ImageFinishGadget.cpp
								 CropAndCombineGadget.cpp
								 ImageWriterGadget.cpp
								 MRIAcquisitionReader.cpp
								 MRIImageWriter.cpp
								 NoiseAdjustGadget.cpp
								 ExtractGadget.cpp
								 FloatToUShortGadget.cpp
								 RemoveROOversamplingGadget.cpp
								 CoilReductionGadget.cpp
								 AutoScaleGadget.cpp)

target_link_libraries(gadgetroncore tinyxml cpucore ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (FILES	AccumulatorGadget.h
			AcquisitionFinishGadget.h
			AcquisitionPassthroughGadget.h
			CropAndCombineGadget.h	
			ExtractGadget.h
			FloatToUShortGadget.h
			FFTGadget.h
			GadgetMRIHeaders.h
			ImageFinishGadget.h
			ImageWriterGadget.h
			MRIAcquisitionReader.h
			MRIImageWriter.h
			NoiseAdjustGadget.h
			RemoveROOversamplingGadget.h
			CoilReductionGadget.h
			AutoScaleGadget.h
			gadgetroncore_export.h
			DESTINATION include)

install(TARGETS gadgetroncore DESTINATION lib)
install(FILES default.xml default_short.xml DESTINATION config)
