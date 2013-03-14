include_directories(${GTEST_INCLUDE})

include_directories(${CMAKE_SOURCE_DIR}/toolboxes/cpucore
					${CMAKE_SOURCE_DIR}/toolboxes/gpucore
					${CUDA_INCLUDE_DIRS})

cuda_add_executable(test_all tests.cpp cuGTBLAS_test.cpp cuNDArray_test.cpp hoCuGTBLAS_test.cpp hoGTBLAS_test.cpp hoNDArray_test hoNDArray_utils_test.cpp 
cuNDArray_Vector_td_test.cpp 
 cuVector_td_test_kernels.cu vector_td_test.cpp)
target_link_libraries(test_all cpucore gpucore ${GTEST_LIBRARIES} gpunfft hostutils ${CUDA_LIBRARIES} pthread )

add_test(test_all test_all)
