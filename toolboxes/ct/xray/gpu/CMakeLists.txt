find_package(HDF5 REQUIRED HL)

  if(${CUDA_VERSION} VERSION_GREATER "4.99")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3} ${CUDA_NVCC_FLAGS4})
  else(${CUDA_VERSION} VERSION_GREATER "4.99")    

    if(${CUDA_VERSION} VERSION_GREATER "4.1")
      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3})
    else(${CUDA_VERSION} VERSION_GREATER "4.1")      

      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2})          

    endif(${CUDA_VERSION} VERSION_GREATER "4.1")
  endif(${CUDA_VERSION} VERSION_GREATER "4.99")

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

cuda_add_library(gpuxray SHARED
  conebeam_projection.cu 
  hoCuConebeamProjectionOperator.cpp 
  )

target_link_libraries(
  gpuxray gpucore gpunfft
  ${CUDA_LIBRARIES}
)

install(TARGETS gpuxray DESTINATION lib)

install(FILES 
  CBCT_acquisition.h
  CBCT_binning.h
  conebeam_projection.h
  hoCuConebeamProjectionOperator.h
  gpuxray_export.h 
  DESTINATION include)
