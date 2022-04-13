add_executable(register_CK_3d_gpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_gpu
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_gpu 
  ${CUDA_LIBRARIES}
  )

install(TARGETS register_CK_3d_gpu DESTINATION bin COMPONENT main)
