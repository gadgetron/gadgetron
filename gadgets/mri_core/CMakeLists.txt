IF (WIN32)
	ADD_DEFINITIONS(-DTIXML_USE_STL)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_MRICORE__)
	ADD_DEFINITIONS(-D__BUILD_GADGETS__)
ENDIF (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   ${ACE_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR}
				       ${FFTW3_INCLUDE_DIR}
				       ${ISMRMRD_INCLUDE_DIR}
					   ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/linalg
				       ${CMAKE_SOURCE_DIR}/gadgets/core
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

add_library(gadgetronmricore SHARED MRINoiseAdjustGadget.cpp PCACoilGadget.cpp)

target_link_libraries(gadgetronmricore tinyxml cpucore linalg ${ISMRMRD_LIBRARIES} ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (FILES MRINoiseAdjustGadget.h PCACoilGadget.h DESTINATION include)

install(TARGETS gadgetronmricore DESTINATION lib)
install(FILES default_optimized.xml radial_single_optimized.xml radial_ms_optimized.xml DESTINATION config)
