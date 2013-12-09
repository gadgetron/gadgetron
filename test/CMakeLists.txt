ENABLE_TESTING()

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${Boost_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  ${ARMADILLO_INCLUDE_DIRS}
  ${GTEST_INCLUDE_DIRS}
  )

if ( CUDA_FOUND )
    cuda_add_executable(test_all 
      tests.cpp 
      hoNDArray_elemwise_test.cpp 
      hoNDArray_operators_test.cpp 
      hoNDArray_blas_test.cpp 
      hoNDArray_utils_test.cpp
      vector_td_test.cpp
      cuNDArray_elemwise_test.cpp 
      cuNDArray_operators_test.cpp 
      cuNDArray_blas_test.cpp 
      cuNDArray_utils_test.cpp
      vector_td_test.cpp
      cuVector_td_test_kernels.h 
      cuVector_td_test_kernels.cu 
      )
else ( CUDA_FOUND )
    add_executable(test_all 
      tests.cpp 
      hoNDArray_elemwise_test.cpp 
      hoNDArray_operators_test.cpp 
      hoNDArray_blas_test.cpp 
      hoNDArray_utils_test.cpp
      )
endif ( CUDA_FOUND )

target_link_libraries(test_all 
  cpucore 
  cpucore_math
  ${BOOST_LIBRARIES}
  ${GTEST_LIBRARIES} 
  )

if ( CUDA_FOUND )

    target_link_libraries(test_all 
      gpucore
      ${CUDA_LIBRARIES} 
      )

endif ( CUDA_FOUND )

add_test(test_all test_all)
