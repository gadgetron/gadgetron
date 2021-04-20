add_executable(register_HS_2d_gpu register_HS_2d.cpp)
add_executable(register_CGHS_2d_gpu register_CGHS_2d.cpp)
add_executable(register_CK_2d_gpu register_CK_2d.cpp)
#add_executable(test_reg_sense_recon test_reg_sense_recon.cpp)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
)

target_link_libraries(register_HS_2d_gpu 
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_gpureg 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpuoperators
  gadgetron_toolbox_gpusolvers 
  ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES}
  )

target_link_libraries(register_CK_2d_gpu
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_gpureg 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpuoperators
  gadgetron_toolbox_gpusolvers 
  ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES}
  )

#target_link_libraries(test_reg_sense_recon 
#  hostutils 
#  gpureg 
#  gpucore 
#  gpuoperators
#  gpusolvers 
#  gpunfft 
#  gpuparallelmri 
#  ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES}
#  )
  
  target_link_libraries(register_CGHS_2d_gpu 
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_gpureg 
  gadgetron_toolbox_gpucore 
  gadgetron_toolbox_gpuoperators
  gadgetron_toolbox_gpusolvers 
  ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES}
  )

install(TARGETS 
  register_HS_2d_gpu
  register_CGHS_2d_gpu  
  register_CK_2d_gpu 
  DESTINATION bin COMPONENT main)
