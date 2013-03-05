include_directories(			${CUDA_INCLUDE_DIRS}
				       ${Boost_INCLUDE_DIR} 
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/operators
						${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
						${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

add_executable(radial_ktsense  main.cpp)

target_link_libraries(radial_ktsense gpucore gpuparallelmri operators gpunfft hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS radial_ktsense DESTINATION bin)
