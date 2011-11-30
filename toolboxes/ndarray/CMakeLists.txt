include_directories(${ACE_INCLUDE_DIR} 
				    ${Boost_INCLUDE_DIR}
				    ${FFTW3_INCLUDE_DIR}  
			        ${CMAKE_SOURCE_DIR}/apps/gadgetron)


add_library(hondarray SHARED FFT.cpp)

target_link_libraries(hondarray ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY})

install(TARGETS hondarray DESTINATION lib)

install(FILES 	
	NDArray.h
	hoNDArray.h
	hoNDArray_fileio.h
	ndarray_export.h
	FFT.h
DESTINATION include)
