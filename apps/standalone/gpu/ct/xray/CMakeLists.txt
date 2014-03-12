find_package(HDF5 REQUIRED HL)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

add_executable(CBCT_reconstruct_FDK_3d CBCT_reconstruct_FDK_3d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_3d gpuxray cpucore cpucore_math cpureg gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_FDK_4d CBCT_reconstruct_FDK_4d.cpp)
target_link_libraries(CBCT_reconstruct_FDK_4d gpuxraycpucore cpucore_math cpureg gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_CG CBCT_reconstruct_CG.cpp)
target_link_libraries(CBCT_reconstruct_CG gpuxray cpucore cpucore_math gpucore cpureg gpuoperators hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_NLCG CBCT_reconstruct_NLCG.cpp)
target_link_libraries(CBCT_reconstruct_NLCG gpuxraycpucore cpucore_math gpucore cpureg gpuoperators hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_reconstruct_SB CBCT_reconstruct_SB.cpp)
target_link_libraries(CBCT_reconstruct_SB gpuxray cpucore cpucore_math gpucore cpureg gpuoperators hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_forwards_projection CBCT_forwards_projection.cpp)
target_link_libraries(CBCT_forwards_projection gpuxraycpucore cpucore_math gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

install(TARGETS 
  CBCT_reconstruct_FDK_3d
  CBCT_reconstruct_FDK_4d
  CBCT_reconstruct_CG
  CBCT_reconstruct_NLCG
  CBCT_reconstruct_SB
  CBCT_forwards_projection
  DESTINATION bin)
