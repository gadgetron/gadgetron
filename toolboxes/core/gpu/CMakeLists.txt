if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUCORE__)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()


cuda_add_library(gadgetron_toolbox_gpucore SHARED
    check_CUDA.h
    CUBLASContextProvider.h
    cudaDeviceManager.h
    cuNDArray.h
    cuNDArray_blas.h
    cuNDArray_elemwise.h
    cuNDArray_operators.h
    cuNDArray_utils.h
    cuNDArray_fileio.h
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
    cuSparseMatrix.cu
  )

set_target_properties(gadgetron_toolbox_gpucore PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpucore 
  gadgetron_toolbox_cpucore
  gadgetron_toolbox_log
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES} 
  ${CUDA_CUFFT_LIBRARIES} 
  ${CUDA_CUBLAS_LIBRARIES} 
  ${CUDA_CUSPARSE_LIBRARIES}
  ${MKL_LIBRARIES}
  )

target_include_directories(gadgetron_toolbox_gpucore
        PUBLIC
		$<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
 

install(TARGETS gadgetron_toolbox_gpucore
  EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
  gpucore_export.h
  cuNDArray.h
  cuNDArray_operators.h
  cuNDArray_elemwise.h
  cuNDArray_blas.h
  cuNDArray_utils.h
  cuNDArray_math.h
  cuNDArray_fileio.h
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
  cuSparseMatrix.h
  cuNDArray_converter.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

  set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)