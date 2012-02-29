if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUNFFT__)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)


include_directories(			${CUDA_INCLUDE_DIRS}
					${Boost_INCLUDE_DIR}
				       ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				       ${CMAKE_SOURCE_DIR}/toolboxes/gpucore)

cuda_add_library(gpunfft SHARED NFFT.cu)

target_link_libraries(gpunfft gpucore ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

install(TARGETS gpunfft DESTINATION lib)

install(FILES NFFT.h gpunfft_export.h DESTINATION include)


