if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

add_executable(blur_3d blur_3d.cpp)
add_executable(deblur_3d_cg deblur_3d_cg.cpp)
add_executable(deblur_3d_sb deblur_3d_sb.cpp)

target_link_libraries(deblur_3d_cg gpucore hostutils gpuoperators gpusolvers gpunfft ${CUDA_LIBRARIES})
target_link_libraries(deblur_3d_sb gpucore hostutils gpuoperators gpusolvers gpunfft ${CUDA_LIBRARIES})
target_link_libraries(blur_3d gpucore hostutils gpuoperators gpusolvers gpunfft ${CUDA_LIBRARIES})

install(TARGETS blur_3d deblur_3d_cg deblur_3d_sb DESTINATION bin)
