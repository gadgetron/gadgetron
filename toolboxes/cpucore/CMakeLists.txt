if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
		    ${Boost_INCLUDE_DIR}
		    ${FFTW3_INCLUDE_DIR})

configure_file( device_config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/device_config.h )

add_library(cpucore SHARED FFT.cpp)

target_link_libraries(cpucore ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS cpucore DESTINATION lib)

install(FILES 	
	NDArray.h
	hoNDArray.h
	hoNDArray_fileio.h
	cpucore_export.h
	FFT.h
	GadgetronTimer.h
	complext.h
	vector_td.h
	vector_td_operators.h
	vector_td_utilities.h
	real_utilities.h
	device_config.h
	DESTINATION include)
