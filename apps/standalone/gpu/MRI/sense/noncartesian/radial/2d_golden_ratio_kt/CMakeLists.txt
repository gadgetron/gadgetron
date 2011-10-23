include_directories(			${CUDA_INCLUDE_DIRS}
						${CUDA_SDK_ROOT_DIR}/common/inc
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/gpuct
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(radial_ktsense  main.cpp)

target_link_libraries(radial_ktsense gpucore gpuparallelmri gpunfft hostutils gpuct solvers ${CUDA_LIBRARIES})

install(TARGETS radial_ktsense DESTINATION bin)
