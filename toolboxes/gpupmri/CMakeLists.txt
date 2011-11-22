if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUPMRI__)
ADD_DEFINITIONS(-DWIN32_LEAN_AND_MEAN)
endif (WIN32)

find_package(CULA REQUIRED)

include_directories(		${CUDA_INCLUDE_DIRS}
				${CULA_INCLUDE_DIR} 
				       ${Boost_INCLUDE_DIR} 
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

target_link_libraries(gpuparallelmri gpucore gpunfft ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CULA_LIBRARIES})

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
	gpupmri_export.h
DESTINATION include)

