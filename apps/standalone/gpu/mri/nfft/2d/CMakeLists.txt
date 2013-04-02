add_executable(nfft main_nfft.cpp)
add_executable(nffth main_nffth.cpp)
add_executable(nffth_cg main_cg.cpp)
add_executable(nffth_sb main_sb.cpp)

target_link_libraries(nfft gpucore gpuoperators gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth gpucore gpuoperators  gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_cg gpusolvers gpuoperators gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_sb gpusolvers gpuoperators gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft nffth nffth_cg nffth_sb DESTINATION bin)
