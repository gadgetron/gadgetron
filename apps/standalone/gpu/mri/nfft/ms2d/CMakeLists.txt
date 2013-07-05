add_executable(nfft_ms nfft_main.cpp)
add_executable(nffth_ms nffth_main.cpp)

target_link_libraries(nfft_ms gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_ms gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft_ms nffth_ms DESTINATION bin)
