if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories(	${CUDA_INCLUDE_DIRS}
					    ${Boost_INCLUDE_DIR}
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(denoise_L2 denoise_L2.cpp)
add_executable(denoise_TV denoise_TV.cpp)

target_link_libraries(denoise_L2 gpucore hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(denoise_TV gpucore hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS denoise_L2 denoise_TV DESTINATION bin)
