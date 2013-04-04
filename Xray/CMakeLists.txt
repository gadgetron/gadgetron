find_package(HDF5 REQUIRED HL)


include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${GADGETRON_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/solvers
  ${CMAKE_SOURCE_DIR}/operators
)

cuda_add_library(Xray_kernels SHARED
  conebeam_projection.cu hoCudaConebeamProjectionOperator.cpp
)

target_link_libraries(Xray_kernels gpucore gpunfft hostutils ${CUDA_LIBRARIES})

add_executable(FDK_reconstruct FDK_reconstruct.cpp )
target_link_libraries(FDK_reconstruct Xray_kernels cpucore cpucore_math gpucore hostutils  ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})
add_executable(ProjectImage ProjectImage.cpp )
target_link_libraries(ProjectImage Xray_kernels cpucore cpucore_math gpucore hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})
#install(TARGETS Reconstruct DESTINATION bin)
