
find_package(CULA REQUIRED)

include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
					${CULA_INCLUDE_DIR} 
				       ${TINYXML_INCLUDE_DIR} 
				       ${CMAKE_SOURCE_DIR}/apps/gadgetron 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
				       ${CMAKE_SOURCE_DIR}/toolboxes/solvers
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
				       ${CMAKE_SOURCE_DIR}/gadgets/core)

cuda_add_library(gpuparallelmri SHARED 	b1_map.cu
									cuSenseOperator.cu
									cuCartesianSenseOperator.cu
									cuNonCartesianSenseOperator.cu
									cuNonCartesianKtSenseOperator.cu
									cuSenseRHSBuffer.cu
									htgrappa.cu)

target_link_libraries(gpuparallelmri gpucore gpunfft ${ACE_LIBRARIES}
${TINYXML_LIBRARIES} ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CULA_LIBRARIES})

install(TARGETS gpuparallelmri DESTINATION lib)

install(FILES 
	b1_map.h
	htgrappa.h
	senseOperator.h
	cuSenseOperator.h
	cuCartesianSenseOperator.h
	cuNonCartesianSenseOperator.h
	cuNonCartesianKtSenseOperator.h
	cuSenseRHSBuffer.h
DESTINATION include)

