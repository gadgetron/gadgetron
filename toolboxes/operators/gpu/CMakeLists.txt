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
  cuImageOperator.h
  cuDiagonalOperator.h
  cuPartialDerivativeOperator.h
  cuConvolutionOperator.h
  cuLaplaceOperator.h
  cuTvOperator.h
  cuTvPicsOperator.h
  cuTv1dOperator.h
  cuDownsampleOperator.h
  cuUpsampleOperator.h
  hoCuTvOperator.h
  hoCuTvPicsOperator.h
  gpuoperators_export.h
  hoCuOperator.h
  DESTINATION include)
