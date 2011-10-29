include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

cuda_add_library(gpunfft SHARED NFFT.cu)

target_link_libraries(gpunfft gpucore ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

install(TARGETS gpunfft DESTINATION lib)

install(FILES NFFT.h DESTINATION include)


