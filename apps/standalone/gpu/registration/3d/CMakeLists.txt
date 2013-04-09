add_executable(register_CK_3d register_CK_3d.cpp)

target_link_libraries(register_CK_3d 
  hostutils 
  gpureg 
  gpucore 
  gpuoperators
  gpusolvers 
  ${CUDA_LIBRARIES}
  )

install(TARGETS register_CK_3d DESTINATION bin)
