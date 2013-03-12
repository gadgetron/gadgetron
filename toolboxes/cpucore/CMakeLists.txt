if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUCORE__)
endif (WIN32)

include_directories(${ACE_INCLUDE_DIR} 
		    ${Boost_INCLUDE_DIR}
		    ${FFTW3_INCLUDE_DIR})

add_library(cpucore SHARED FFT.cpp hoGTBLAS.cpp)

target_link_libraries(cpucore ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ARMADILLO_LIBRARIES})

install(TARGETS cpucore DESTINATION lib)

install(FILES 	
	NDArray.h
	hoNDArray.h
	hoNDArray_fileio.h
	GadgetronException.h
	cpucore_export.h
	FFT.h
	GadgetronTimer.h
	complext.h
	vector_td.h
	vector_td_operators.h
	vector_td_utilities.h
	real_utilities.h
	cpucore_defines.h
	hoGTBLAS.h
	hoNDArray_utils.h
	gtArmadillo.h
	DESTINATION include)
