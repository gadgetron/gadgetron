if (GTEST_FOUND AND ARMADILLO_FOUND)

ENABLE_TESTING()

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/wavelet/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${Boost_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
  ${GTEST_INCLUDE_DIRS}
  ${ACE_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  )

if (CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
endif (CUDA_FOUND)

link_libraries(
    gadgetron_toolbox_cpucore 
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpuwavelet
    ${BOOST_LIBRARIES}
    ${GTEST_LIBRARIES} 
    ${ARMADILLO_LIBRARIES}
    )

if ( CUDA_FOUND )

    include_directories( ${CUDA_INCLUDE_DIRS} )

  cuda_add_executable(test_all 
      tests.cpp 
      hoNDArray_elemwise_test.cpp 
      hoNDArray_blas_test.cpp 
      hoNDArray_utils_test.cpp
      hoNDFFT_test.cpp
      hoNDWavelet_test.cpp
      vector_td_test.cpp
      cuNDArray_elemwise_test.cpp 
      cuNDArray_operators_test.cpp 
      cuNDArray_blas_test.cpp 
      cuNDArray_utils_test.cpp
      vector_td_test.cpp
      cuVector_td_test_kernels.h 
      cuVector_td_test_kernels.cu 
      cuNDFFT_test.cpp
      )
else ( CUDA_FOUND )
    add_executable(test_all 
      tests.cpp 
      hoNDArray_elemwise_test.cpp 
      hoNDArray_blas_test.cpp 
      hoNDArray_utils_test.cpp
      hoNDFFT_test.cpp
      hoNDWavelet_test.cpp
      )
endif ( CUDA_FOUND )

if ( CUDA_FOUND )
  target_link_libraries(test_all 
    gadgetron_toolbox_gpucore
    gadgetron_toolbox_gpufft
    )
endif( CUDA_FOUND )

add_test(test_all test_all)

endif (GTEST_FOUND AND ARMADILLO_FOUND)

add_subdirectory(integration)
