if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

add_executable(blur_2d blur_2d.cpp)
add_executable(deblur_2d_cg deblur_2d_cg.cpp)
add_executable(deblur_2d_sb deblur_2d_sb.cpp)

target_link_libraries(deblur_2d_cg gpucore hostutils gpuoperators gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(deblur_2d_sb gpucore hostutils gpuoperators gpusolvers ${CUDA_LIBRARIES})
target_link_libraries(blur_2d gpucore hostutils gpuoperators gpusolvers ${CUDA_LIBRARIES})

install(TARGETS blur_2d deblur_2d_cg deblur_2d_sb DESTINATION bin)

