if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUCORE__)
endif (WIN32)

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( 
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
)

cuda_add_library(gadgetron_toolbox_gpucore SHARED 
    check_CUDA.h
    CUBLASContextProvider.h
    cudaDeviceManager.h
    cuNDArray.h
    cuNDArray_blas.h
    cuNDArray_elemwise.h
    cuNDArray_operators.h
    cuNDArray_utils.h
    cuNDArray_reductions.h
    GadgetronCuException.h
    gpucore_export.h
    GPUTimer.h
    hoCuNDArray.h
    hoCuNDArray_blas.h
    hoCuNDArray_elemwise.h
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
    hoCuNDArray_blas.cpp
    CUBLASContextProvider.cpp
    cudaDeviceManager.cpp
  )

set_target_properties(gadgetron_toolbox_gpucore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpucore 
  gadgetron_toolbox_cpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES} 
  ${MKL_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_gpucore DESTINATION lib COMPONENT main)

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
  hoCuNDArray_elemwise.h
  hoCuNDArray_utils.h
  hoCuNDArray_math.h
  GPUTimer.h				
  GadgetronCuException.h
  radial_utilities.h
  real_utilities_device.h
  check_CUDA.h
  cudaDeviceManager.h
  CUBLASContextProvider.h
  setup_grid.h
  DESTINATION include COMPONENT main)
