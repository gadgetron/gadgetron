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

add_executable(blur_3d blur_3d.cpp)
add_executable(deblur_3d_cg deblur_3d_cg.cpp)
add_executable(deblur_3d_sb deblur_3d_sb.cpp)

target_link_libraries(deblur_3d_cg gpucore hostutils operators solvers ${CUDA_LIBRARIES})
target_link_libraries(deblur_3d_sb gpucore hostutils operators solvers ${CUDA_LIBRARIES})
target_link_libraries(blur_3d gpucore hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS blur_3d deblur_3d_cg deblur_3d_sb DESTINATION bin)
