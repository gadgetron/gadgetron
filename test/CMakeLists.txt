include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${GTEST_INCLUDE_DIRS}
  )

if(CUDA_FOUND)
  include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CUDA_INCLUDE_DIRS}
    )
endif(CUDA_FOUND)

cuda_add_executable(test_all 
  tests.cpp 
  hoNDArray_elemwise_test.cpp 
  hoNDArray_operators_test.cpp 
  hoNDArray_blas_test.cpp 
  hoNDArray_utils_test.cpp
  cuNDArray_elemwise_test.cpp 
  cuNDArray_operators_test.cpp 
  cuNDArray_blas_test.cpp 
  cuNDArray_utils_test.cpp
#  hoCuGTBLAS_test.cpp 
#  cuNDArray_Vector_td_test.cpp 
  vector_td_test.cpp
  cuVector_td_test_kernels.cu 
  )

target_link_libraries(test_all 
  cpucore 
  cpucore_math
  gpucore
  pthread
  ${BOOST_LIBRARIES}	
  ${GTEST_LIBRARIES} 
  )

if(CUDA_FOUND)
  target_link_libraries(test_all 
    #gpucore 
    #gpunfft 
    ${CUDA_LIBRARIES} 
    )
endif(CUDA_FOUND)

add_test(test_all test_all)
