
include_directories(		${ACE_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
					${CUDA_SDK_ROOT_DIR}/common/inc 
				       ${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_library(cutil STATIC IMPORTED)
set_target_properties(cutil PROPERTIES IMPORTED_LOCATION ${CUDA_SDK_ROOT_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}cutil_x86_64${CMAKE_STATIC_LIBRARY_SUFFIX})

cuda_add_library(gpucore SHARED ndarray_vector_td_utilities.cu
							cuNDArray_kernels.cu
							radial_utilities.cu
							cuNDArray.cpp
							cuNDFFT.cpp)

target_link_libraries(gpucore ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES}
${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} cutil)

install(TARGETS gpucore DESTINATION lib)

install(FILES
GPUTimer.h				
check_CUDA.h
vector_td.h
vector_td_operators.h
vector_td_utilities.h
ndarray_vector_td_utilities.h
cuNDArray.h
cuNDFFT.h
radial_utilities.h
real_utilities.h
real_utilities_device.h
DESTINATION include)
