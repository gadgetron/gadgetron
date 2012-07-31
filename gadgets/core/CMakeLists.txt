IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CORE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

find_path(GADGETRON_CMAKE_MODULES FindGadgetron.cmake HINTS
$ENV{GADGETRON_HOME}/cmake
/usr/local/gadgetron/cmake)

if (NOT GADGETRON_CMAKE_MODULES)
   #TODO: Check if path found otherwise create error
   MESSAGE(FATAL_ERROR "GADGETRON_CMAKE_MODULES cannot be found. 
   Try to set GADGETRON_HOME environment variable.")
endif(NOT GADGETRON_CMAKE_MODULES)

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES})

find_package(Ismrmrd REQUIRED)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
				       ${ISMRMRD_INCLUDE_DIR}
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
								 AutoScaleGadget.cpp
								 GadgetIsmrmrdReadWrite.cpp)

target_link_libraries(gadgetroncore tinyxml cpucore ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

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
			GadgetIsmrmrdReadWrite.h
			gadgetroncore_export.h
			DESTINATION include)

install(TARGETS gadgetroncore DESTINATION lib)
install(FILES default.xml default_short.xml DESTINATION config)
