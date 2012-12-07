
if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUCORE__)
endif (WIN32)

find_package(CULA REQUIRED)

include_directories( ${CUDA_INCLUDE_DIRS}
                     ${CULA_INCLUDE_DIR}
                     ${Boost_INCLUDE_DIR}
                     ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
                     ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
                     ${CMAKE_SOURCE_DIR}/toolboxes/hostutils )

add_executable(cufft_test cufft_test.cpp)

cuda_add_library(gpucore SHARED 
                            ndarray_vector_td_utilities.cu
                            cuNDArray_kernels.cu
                            radial_utilities.cu
                            cuNDArray.cu
                            cuNDFFT.cpp
                            cuGTBLAS.cu
                            CUBLASContextProvider.cpp
							)

target_link_libraries(gpucore ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CULA_LIBRARIES})
target_link_libraries(cufft_test ${FFTW3_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CULA_LIBRARIES} gpucore gpucore)

install(TARGETS gpucore DESTINATION lib)
install(TARGETS cufft_test DESTINATION bin)

install(FILES
GPUTimer.h				
check_CUDA.h
ndarray_vector_td_utilities.h
cuNDArray.h
cuNDFFT.h
cuGTBLAS.h
GadgetronCuException.h
radial_utilities.h
real_utilities_device.h
gpucore_export.h
CUBLASContextProvider.h
DESTINATION include)
