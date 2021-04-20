add_executable(register_CK_3d_gpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_gpu
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_gpureg 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpuoperators
  gadgetron_toolbox_gpusolvers 
  ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES}
  )

install(TARGETS register_CK_3d_gpu DESTINATION bin COMPONENT main)
