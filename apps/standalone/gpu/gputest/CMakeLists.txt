
include_directories(			${CUDA_INCLUDE_DIRS}
				       ${Boost_INCLUDE_DIR} 
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/gpupmri 
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(gputestapp  main.cpp)

target_link_libraries(gputestapp gpucore gpunfft gpuparallelmri hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS gputestapp DESTINATION bin)
