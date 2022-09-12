if (WIN32)
add_definitions(-D_USE_MATH_DEFINES)
endif ()

add_executable(blur_3d blur_3d.cpp)
add_executable(deblur_3d_cg deblur_3d_cg.cpp)
add_executable(deblur_3d_sb deblur_3d_sb.cpp)

target_link_libraries(deblur_3d_cg gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(deblur_3d_sb gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})
target_link_libraries(blur_3d gadgetron_toolbox_gpu gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})

install(TARGETS blur_3d deblur_3d_cg deblur_3d_sb DESTINATION bin COMPONENT main)
