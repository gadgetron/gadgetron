add_executable(nfft_2d_ms_radial nfft_main.cpp)
add_executable(nffth_2d_ms_radial nffth_main.cpp)
add_executable(nffth_2d_ms_generic nffth_generic.cpp)

target_link_libraries(nfft_2d_ms_radial gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_2d_ms_radial gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_2d_ms_generic gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft_2d_ms_radial nffth_2d_ms_radial nffth_2d_ms_generic DESTINATION bin COMPONENT main)
