if (WIN32)
add_definitions(-D_USE_MATH_DEFINES)
endif ()

add_executable(blur_2d blur_2d.cpp)
add_executable(deblur_2d_cg deblur_2d_cg.cpp)
add_executable(deblur_2d_sb deblur_2d_sb.cpp)

target_link_libraries(deblur_2d_cg gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(deblur_2d_sb gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(blur_2d gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})

install(TARGETS blur_2d deblur_2d_cg deblur_2d_sb DESTINATION bin COMPONENT main)

