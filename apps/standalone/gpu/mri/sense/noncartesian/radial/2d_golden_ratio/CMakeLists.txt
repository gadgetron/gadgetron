add_executable(sense_cg_radial_2d main_cg.cpp)
add_executable(sense_sbc_radial_2d main_sbc.cpp)
add_executable(sense_gpbb_radial_2d main_gpbb.cpp)
add_executable(sense_nlcg_radial_2d main_nlcg.cpp)

target_link_libraries(sense_cg_radial_2d gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(sense_sbc_radial_2d gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(sense_gpbb_radial_2d gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(sense_nlcg_radial_2d gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

install(TARGETS sense_cg_radial_2d sense_sbc_radial_2d sense_gpbb_radial_2d sense_nlcg_radial_2d DESTINATION bin)
