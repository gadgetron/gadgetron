add_executable(radial_sense_cg main_cg.cpp)
add_executable(radial_sense_sbc main_sbc.cpp)

target_link_libraries(radial_sense_cg gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(radial_sense_sbc gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

install(TARGETS radial_sense_cg radial_sense_sbc DESTINATION bin)
