find_package(HDF5 REQUIRED HL)

  if(${CUDA_VERSION} VERSION_GREATER "4.99")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3} ${CUDA_NVCC_FLAGS4})
  else()    

    if(${CUDA_VERSION} VERSION_GREATER "4.1")
      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2} ${CUDA_NVCC_FLAGS3})
    else()

      set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS2})

    endif()
  endif()

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/ct/xray/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/gpu
    ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
    ${CUDA_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}/cpp
    ${ISMRMRD_INCLUDE_DIR}
)

link_directories(${Boost_LIBRARY_DIR} ${HDF5_LIBRARY_DIRS})
link_libraries(${HDF5_LIBRARIES})

add_executable(CBCT_reconstruct_FDK_3d CBCT_reconstruct_FDK_3d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_3d gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_reconstruct_FDK_4d CBCT_reconstruct_FDK_4d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_4d gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_reconstruct_CG CBCT_reconstruct_CG.cpp)
target_link_libraries(CBCT_reconstruct_CG gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_gpuoperators gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_reconstruct_NLCG CBCT_reconstruct_NLCG.cpp)
target_link_libraries(CBCT_reconstruct_NLCG gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_gpuoperators gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_reconstruct_SB CBCT_reconstruct_SB.cpp)
target_link_libraries(CBCT_reconstruct_SB gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_gpuoperators gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_forwards_projection CBCT_forwards_projection.cpp)
target_link_libraries(CBCT_forwards_projection gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})

add_executable(CBCT_OF_reconstruct CBCT_OF_reconstruct.cpp)
target_link_libraries(CBCT_OF_reconstruct gadgetron_toolbox_gpuxray gadgetron_toolbox_cpucore gadgetron_toolbox_cpucore_math gadgetron_toolbox_gpureg gadgetron_toolbox_cpureg gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})


install(TARGETS 
  CBCT_reconstruct_FDK_3d
  CBCT_reconstruct_FDK_4d
  CBCT_reconstruct_CG
  CBCT_reconstruct_NLCG
  CBCT_reconstruct_SB
  CBCT_forwards_projection
  DESTINATION bin COMPONENT main)
