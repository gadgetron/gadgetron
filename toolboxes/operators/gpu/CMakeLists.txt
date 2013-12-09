if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUOPERATORS__)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

cuda_add_library(gpuoperators SHARED 
    cuConvolutionOperator.h
    cuDiagonalOperator.h
    cuFFTOperator.h
    cuIdentityOperator.h
    cuImageOperator.h
    cuLaplaceOperator.h
    cuMultiplicationOperatorContainer.h
    cuPartialDerivativeOperator.h
    cuTv1dOperator.h
    cuTvOperator.h
    cuTvPicsOperator.h
    gpuoperators_export.h
    hoCuOperator.h
    cuPartialDerivativeOperator.cu
    cuLaplaceOperator.cu
    cuTvOperator.cu
    cuTv1dOperator.cu
    cuConvolutionOperator.cu
  )

target_link_libraries(gpuoperators 
  gpucore 
  gpunfft
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gpuoperators DESTINATION lib)

install(FILES 	
    cuConvolutionOperator.h
    cuDiagonalOperator.h
    cuFFTOperator.h
    cuIdentityOperator.h
    cuImageOperator.h
    cuLaplaceOperator.h
    cuMultiplicationOperatorContainer.h
    cuPartialDerivativeOperator.h
    cuTv1dOperator.h
    cuTvOperator.h
    cuTvPicsOperator.h
    hoCuOperator.h
  DESTINATION include)
