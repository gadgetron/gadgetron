
include_directories( ${CUDA_INCLUDE_DIRS}
		     ${Boost_INCLUDE_DIR} 
		     ${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
		     ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
		     ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
		     ${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

add_executable(nfft nfft_main.cpp)
add_executable(nffth nffth_main.cpp)

target_link_libraries(nfft gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft nffth DESTINATION bin)
