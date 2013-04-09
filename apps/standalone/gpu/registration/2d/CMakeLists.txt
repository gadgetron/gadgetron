add_executable(register_HS_2d register_HS_2d.cpp)
add_executable(register_CK_2d register_CK_2d.cpp)
add_executable(test_reg_sense_recon test_reg_sense_recon.cpp)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
)

target_link_libraries(register_HS_2d 
  hostutils 
  gpureg 
  gpucore 
  gpuoperators
  gpusolvers 
  ${CUDA_LIBRARIES}
  )

target_link_libraries(register_CK_2d 
  hostutils 
  gpureg 
  gpucore 
  gpuoperators
  gpusolvers 
  ${CUDA_LIBRARIES}
  )

target_link_libraries(test_reg_sense_recon 
  hostutils 
  gpureg 
  gpucore 
  gpuoperators
  gpusolvers 
  gpunfft 
  gpuparallelmri 
  ${CUDA_LIBRARIES}
  )

install(TARGETS 
  register_HS_2d 
  register_CK_2d 
  DESTINATION bin)
