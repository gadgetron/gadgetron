include_directories(${GTEST_INCLUDE})

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
  ${CMAKE_SOURCE_DIR}/toolboxes/cpucore/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/gpucore
  ${CUDA_INCLUDE_DIRS}
  )

cuda_add_executable(test_all 
  tests.cpp 
  hoNDArray_elemwise_test.cpp 
  hoNDArray_operators_test.cpp 
  hoNDArray_blas_test.cpp 
#  hoNDArray_utils_test.cpp 
#  cuGTBLAS_test.cpp 
#  cuNDArray_test.cpp 
#  hoCuGTBLAS_test.cpp 
#  cuNDArray_Vector_td_test.cpp 
#  vector_td_test.cpp
#  cuVector_td_test_kernels.cu 
  )

target_link_libraries(test_all 
  cpucore 
  cpucore_math
  #gpucore 
  #gpunfft 
  hostutils 
  ${GTEST_LIBRARIES} 
  ${CUDA_LIBRARIES} 
  pthread
  )

add_test(test_all test_all)
