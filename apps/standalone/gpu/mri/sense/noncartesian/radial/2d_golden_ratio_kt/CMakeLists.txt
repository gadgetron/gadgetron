add_executable(radial_ktsense  main.cpp)

target_link_libraries(radial_ktsense gpucore gpuparallelmri gpuoperators gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

install(TARGETS radial_ktsense DESTINATION bin)
