
include_directories(			${CUDA_INCLUDE_DIRS}
						${CUDA_SDK_ROOT_DIR}/common/inc
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/gpuct
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(nfft  nfft_main.cpp)
add_executable(nffth  nffth_main.cpp)

target_link_libraries(nfft gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft nffth DESTINATION bin)
