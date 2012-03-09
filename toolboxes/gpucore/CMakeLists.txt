
if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUCORE__)
endif (WIN32)

include_directories(		${CUDA_INCLUDE_DIRS}
				${Boost_INCLUDE_DIR}
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore
				${CMAKE_SOURCE_DIR}/toolboxes/gadgettools )

cuda_add_library(gpucore SHARED ndarray_vector_td_utilities.cu
							cuNDArray_kernels.cu
							radial_utilities.cu
							cuNDArray.cpp
							cuNDFFT.cpp)

target_link_libraries(gpucore ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

install(TARGETS gpucore DESTINATION lib)

install(FILES
GPUTimer.h				
check_CUDA.h
complext.h
vector_td.h
vector_td_operators.h
vector_td_utilities.h
ndarray_vector_td_utilities.h
cuNDArray.h
cuNDFFT.h
radial_utilities.h
real_utilities.h
real_utilities_device.h
gpucore_export.h
DESTINATION include)
