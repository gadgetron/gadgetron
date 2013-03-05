if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories(	${CUDA_INCLUDE_DIRS}
					    ${Boost_INCLUDE_DIR}
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/operators
						${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

add_executable(blur_2d blur_2d.cpp)
add_executable(deblur_2d_cg deblur_2d_cg.cpp)
add_executable(deblur_2d_sb deblur_2d_sb.cpp)

target_link_libraries(deblur_2d_cg gpucore hostutils operators solvers ${CUDA_LIBRARIES})
target_link_libraries(deblur_2d_sb gpucore hostutils operators  solvers ${CUDA_LIBRARIES})
target_link_libraries(blur_2d gpucore hostutils operators  solvers ${CUDA_LIBRARIES})

install(TARGETS blur_2d deblur_2d_cg deblur_2d_sb DESTINATION bin)

