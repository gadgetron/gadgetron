IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CORE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

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
INCLUDE_DIRECTORIES(${XSD_INCLUDES} ${XERCESC_INCLUDE_DIR})

add_library(gadgetroncore SHARED AcquisitionPassthroughGadget.cpp
								 AcquisitionFinishGadget.cpp 
								 AccumulatorGadget.cpp
								 FFTGadget.cpp
								 ImageFinishGadget.cpp
								 CropAndCombineGadget.cpp
								 ImageWriterGadget.cpp
								 MRIImageWriter.cpp
								 NoiseAdjustGadget.cpp
								 ExtractGadget.cpp
								 FloatToUShortGadget.cpp
								 RemoveROOversamplingGadget.cpp
								 CoilReductionGadget.cpp
								 AutoScaleGadget.cpp
								 GadgetIsmrmrdReadWrite.cpp
								 ${XSDS_SOURCES}
								 )

target_link_libraries(gadgetroncore 
  tinyxml cpucore 
  ${Boost_LIBRARIES}
  ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} 
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  ${XERCESC_LIBRARIES} 
  )

install (FILES	
            AccumulatorGadget.h
			AcquisitionFinishGadget.h
			AcquisitionPassthroughGadget.h
			CropAndCombineGadget.h	
			ExtractGadget.h
			FloatToUShortGadget.h
			FFTGadget.h
			GadgetMRIHeaders.h
			ImageFinishGadget.h
			ImageWriterGadget.h
			MRIImageWriter.h
			NoiseAdjustGadget.h
			RemoveROOversamplingGadget.h
			CoilReductionGadget.h
			AutoScaleGadget.h
			GadgetIsmrmrdReadWrite.h
			gadgetroncore_export.h
			DESTINATION include)

install(TARGETS gadgetroncore DESTINATION lib)
install(FILES default.xml default_short.xml DESTINATION config)

install(FILES ${ISMRMRD_SCHEMA_DIR}/ismrmrd.xsd DESTINATION schema)
install(FILES ${ISMRMRD_LIBRARIES} DESTINATION lib)

