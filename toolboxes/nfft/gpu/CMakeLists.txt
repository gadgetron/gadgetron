if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUNFFT__)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( 
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

cuda_add_library(gpunfft SHARED 
    cuNFFT.h
    cuNFFTOperator.h
    gpunfft_export.h
    cuNFFT.cu 
    cuNFFTOperator.cu
  )

target_link_libraries(
  gpunfft 
  gpucore 
  ${Boost_LIBRARIES}
  ${FFTW3_LIBRARIES} 
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES}
  )

install(TARGETS gpunfft DESTINATION lib)

install(FILES 
  cuNFFT.h 
  cuNFFTOperator.h 
  gpunfft_export.h 
  DESTINATION include)
