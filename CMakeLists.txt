find_package(HDF5 REQUIRED HL)

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${GADGETRON_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/solvers
  ${CMAKE_SOURCE_DIR}/operators
)

# The Armadillo sparse matrix suppport seems to be too slow be be a real alternative to the eigen version
#

find_path(EIGEN_INCLUDE_DIR Eigen/Sparse 
  HINTS /usr/local/include/eigen3)

if(EIGEN_INCLUDE_DIR)
  set(EIGEN_RESAMPLER hoLinearResampleOperator_eigen.cpp)
  include_directories(${EIGEN_INCLUDE_DIR})  
elseif(EIGEN_INCLUDE_DIR)
  MESSAGE(STATUS "Eigen not found. Not compiling the Eigen based resampler and the CB_OF_reconstruct app") 
endif(EIGEN_INCLUDE_DIR)

cuda_add_library(Xray_kernels SHARED
  conebeam_projection.cu hoCuConebeamProjectionOperator.cpp 
  ${EIGEN_RESAMPLER}
  )

target_link_libraries(Xray_kernels gpucore gpunfft hostutils ${CUDA_LIBRARIES})

add_executable(FDK_reconstruct_3d FDK_reconstruct_3d.cpp)
target_link_libraries(FDK_reconstruct_3d Xray_kernels cpucore cpucore_math cpureg gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(FDK_reconstruct_4d FDK_reconstruct_4d.cpp)
target_link_libraries(FDK_reconstruct_4d Xray_kernels cpucore cpucore_math cpureg gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CB_reconstruct CB_reconstruct.cpp)
target_link_libraries(CB_reconstruct Xray_kernels cpucore cpucore_math gpucore cpureg gpuoperators hostutils  ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CB_SBreconstruct CB_SBreconstruct.cpp)
target_link_libraries(CB_SBreconstruct Xray_kernels cpucore cpucore_math gpucore cpureg gpuoperators hostutils  ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

add_executable(CBCT_forwards_projection CBCT_forwards_projection.cpp)
target_link_libraries(CBCT_forwards_projection Xray_kernels cpucore cpucore_math gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})

if(EIGEN_INCLUDE_DIR)
  add_executable(CB_OF_reconstruct CB_OF_reconstruct.cpp)
  target_link_libraries(CB_OF_reconstruct Xray_kernels cpucore cpucore_math gpucore cpureg gpureg gpuoperators hostutils  ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})
endif(EIGEN_INCLUDE_DIR)

#install(TARGETS Reconstruct DESTINATION bin)
