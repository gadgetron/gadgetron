if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories(	${CUDA_INCLUDE_DIRS}
					    ${Boost_INCLUDE_DIR}
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/ndarray)

add_executable(blur blur.cpp)
add_executable(deblur_cg deblur_cg.cpp)
add_executable(deblur_sb deblur_sb.cpp)

target_link_libraries(deblur_cg gpucore hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(deblur_sb gpucore hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(blur gpucore hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS blur deblur_cg deblur_sb DESTINATION bin)
