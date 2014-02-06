add_executable(ktsense_radial_2d  main.cpp)

target_link_libraries(ktsense_radial_2d gpucore gpuparallelmri gpuoperators gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

install(TARGETS ktsense_radial_2d DESTINATION bin)
