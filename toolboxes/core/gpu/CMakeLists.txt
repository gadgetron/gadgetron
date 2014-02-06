if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUCORE__)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

find_package(CULA REQUIRED)

include_directories( 
  ${CUDA_INCLUDE_DIRS}
  ${CULA_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
)

cuda_add_library(gpucore SHARED 
    check_CUDA.h
    CUBLASContextProvider.h
    cudaDeviceManager.h
    cuNDArray.h
    cuNDArray_blas.h
    cuNDArray_elemwise.h
    cuNDArray_operators.h
    cuNDArray_utils.h
    cuNDArray_reductions.h
    cuNDFFT.h
    cuNDFFT.cpp
    GadgetronCuException.h
    gpucore_export.h
    GPUTimer.h
    hoCuNDArray.h
    hoCuNDArray_blas.h
    hoCuNDArray_elemwise.h
    hoCuNDArray_operators.h
    hoCuNDArray_utils.h
    radial_utilities.h
    real_utilities_device.h
    setup_grid.h
    cuNDArray_operators.cu
    cuNDArray_elemwise.cu
    cuNDArray_blas.cu
    cuNDArray_utils.cu
    cuNDArray_reductions.cu
    radial_utilities.cu
    hoCuNDArray_blas.cu
    CUBLASContextProvider.cpp
    cudaDeviceManager.cpp
  )

target_link_libraries(gpucore cpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES} 
  ${CULA_LIBRARIES}
  )

install(TARGETS gpucore DESTINATION lib)

install(FILES
  gpucore_export.h
  cuNDArray.h
  cuNDArray_operators.h
  cuNDArray_elemwise.h
  cuNDArray_blas.h
  cuNDArray_utils.h
  cuNDArray_math.h
  cuNDArray_reductions.h
  hoCuNDArray.h
  hoCuNDArray_blas.h
  hoCuNDArray_operators.h
  hoCuNDArray_elemwise.h
  hoCuNDArray_utils.h
  hoCuNDArray_math.h
  GPUTimer.h				
  cuNDFFT.h
  GadgetronCuException.h
  radial_utilities.h
  real_utilities_device.h
  check_CUDA.h
  cudaDeviceManager.h
  CUBLASContextProvider.h
  setup_grid.h
  DESTINATION include)
