include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
)

add_executable(nfft_2d_radial main_nfft.cpp)
add_executable(nffth_2d_radial main_nffth.cpp)
#add_executable(moco moco.cpp)
add_executable(nffth_cg_2d_radial main_cg.cpp)
add_executable(nffth_sb_2d_radial main_sb.cpp)
add_executable(test_cuda_gadgetron test_cuda_gadgetron.cpp)

target_link_libraries(nfft_2d_radial gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_2d_radial gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
#target_link_libraries(moco gpusolvers gpureg gpucore gpuparallelmri gpuoperators gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_cg_2d_radial gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_sb_2d_radial gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(test_cuda_gadgetron gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})


install(TARGETS test_cuda_gadgetron nfft_2d_radial nffth_2d_radial nffth_cg_2d_radial nffth_sb_2d_radial DESTINATION bin COMPONENT main)
