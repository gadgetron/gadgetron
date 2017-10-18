if (GTEST_FOUND AND ARMADILLO_FOUND)

enable_testing()

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/image
  ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/image_io
  ${CMAKE_SOURCE_DIR}/toolboxes/cmr
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
  ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/ffd
  ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
  ${CMAKE_SOURCE_DIR}/toolboxes/pattern_recognition
  ${Boost_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
  ${GTEST_INCLUDE_DIRS}
  ${ACE_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  )

if (CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
endif ()

link_libraries(
    gadgetron_toolbox_cpucore 
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpunfft
    gadgetron_toolbox_cpudwt
    gadgetron_toolbox_cpu_image
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_log
    gadgetron_toolbox_cpuklt 
    gadgetron_toolbox_image_analyze_io 
    gadgetron_toolbox_mri_core 
    gadgetron_toolbox_cpuoperator 
    gadgetron_toolbox_cpu_image
    gadgetron_toolbox_cmr
    gadgetron_toolbox_pr
    ${BOOST_LIBRARIES}
    ${GTEST_LIBRARIES} 
    ${ARMADILLO_LIBRARIES}
    )

set(test_src_files tests.cpp 
      hoNDArray_elemwise_test.cpp 
      hoNDArray_blas_test.cpp 
      hoNDArray_utils_test.cpp 
      hoNDArray_reductions_test.cpp 
      hoNDFFT_test.cpp
      hoNFFT_test.cpp
      hoNDWavelet_test.cpp
      curveFitting_test.cpp
      image_morphology_test.cpp 
      pattern_recognition_test.cpp 
      )

if ( CUDA_FOUND )

    include_directories( ${CUDA_INCLUDE_DIRS} )

    cuda_add_executable(test_all 
        ${test_src_files} 
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
else ()
    add_executable(test_all 
        ${test_src_files}
        )
endif ()

if ( CUDA_FOUND )
    target_link_libraries(test_all 
        gadgetron_toolbox_gpucore
        gadgetron_toolbox_gpufft
        )
endif()

add_test(test_all test_all)

endif ()

add_subdirectory(integration)
