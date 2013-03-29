if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUOPERATORS__)
endif (WIN32)

include_directories(		
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

cuda_add_library(gpuoperators SHARED 
  cuPartialDerivativeOperator.cu
  cuLaplaceOperator.cu
  cuTVOperator.cu
  cuTV1DOperator.cu
  cuConvolutionOperator.cu
  )

target_link_libraries(gpuoperators 
  gpucore 
  gpunfft
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gpuoperators DESTINATION lib)

install(FILES 	
  cuImageOperator.h
  cuDiagonalOperator.h
  cuEncodedImageOperator.h
  cuPartialDerivativeOperator.h
  cuConvolutionOperator.h
  cuLaplaceOperator.h
  cuTVOperator.h
  cuTV1DOperator.h
  DESTINATION include)
