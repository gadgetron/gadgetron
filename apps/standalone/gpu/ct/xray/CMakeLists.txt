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
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/ct/xray/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CUDA_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${ISMRMRD_INCLUDE_DIR}
)

link_directories(${Boost_LIBRARY_DIR})
link_libraries(${HDF5_LIBRARIES})

add_executable(CBCT_reconstruct_FDK_3d CBCT_reconstruct_FDK_3d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_3d gpuxray cpucore cpucore_math gpucore hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_FDK_4d CBCT_reconstruct_FDK_4d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_4d gpuxray cpucore cpucore_math gpucore hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_CG CBCT_reconstruct_CG.cpp)
target_link_libraries(CBCT_reconstruct_CG gpuxray cpucore cpucore_math gpucore gpuoperators hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_NLCG CBCT_reconstruct_NLCG.cpp)
target_link_libraries(CBCT_reconstruct_NLCG gpuxray cpucore cpucore_math gpucore gpuoperators hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_SB CBCT_reconstruct_SB.cpp)
target_link_libraries(CBCT_reconstruct_SB gpuxray cpucore cpucore_math gpucore gpuoperators hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_forwards_projection CBCT_forwards_projection.cpp)
target_link_libraries(CBCT_forwards_projection gpuxray cpucore cpucore_math gpucore hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS 
  CBCT_reconstruct_FDK_3d
  CBCT_reconstruct_FDK_4d
  CBCT_reconstruct_CG
  CBCT_reconstruct_NLCG
  CBCT_reconstruct_SB
  CBCT_forwards_projection
  DESTINATION bin COMPONENT main)
