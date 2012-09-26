if (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CGSENSE__)
ADD_DEFINITIONS(-D__BUILD_GADGETS__)
endif (WIN32)

find_package(Ismrmrd REQUIRED)

include_directories(   ${ACE_INCLUDE_DIR}
				       ${CUDA_INCLUDE_DIRS}
				       ${Boost_INCLUDE_DIR}
				       ${ISMRMRD_INCLUDE_DIR}
					   ${XSD_INCLUDE_DIR}
					   ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/solvers
					   ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

add_library(gadgetroncgsense SHARED GPUCGGadget.cpp GPUCGGoldenRadial.cpp GPUCGFixedRadial.cpp GPUCGGadgetGeneric.cpp)

target_link_libraries(gadgetroncgsense solvers gpuparallelmri tinyxml ${ISMRMRD_LIBRARIES} ${ISMRMRD_XSD_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install (TARGETS gadgetroncgsense DESTINATION lib)
install (FILES radial_single.xml  DESTINATION config)
