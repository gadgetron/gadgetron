include_directories(			${CUDA_INCLUDE_DIRS}
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(blur  blur.cpp)
add_executable(cg_deblur  cg_deblur.cpp)
add_executable(sb_deblur  sb_deblur.cpp)

target_link_libraries(cg_deblur gpucore hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(sb_deblur gpucore hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(blur gpucore hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS blur cg_deblur sb_deblur DESTINATION bin)
