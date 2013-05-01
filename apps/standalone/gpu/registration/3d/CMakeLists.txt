add_executable(register_CK_3d_gpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_gpu
  hostutils 
  gpureg 
  gpucore 
  gpuoperators
  gpusolvers 
  ${CUDA_LIBRARIES}
  )

install(TARGETS register_CK_3d_gpu DESTINATION bin)
