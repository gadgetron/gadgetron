if (WIN32)
add_definitions(-D_USE_MATH_DEFINES)
endif ()

include_directories(${CMAKE_SOURCE_DIR}/toolboxes/dwt/gpu)
add_executable(denoise_TV denoise_TV.cpp)

target_link_libraries(denoise_TV gadgetron_toolbox_gpudwt gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils gadgetron_toolbox_gpusolvers gadgetron_toolbox_gpuoperators ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES})

install(TARGETS denoise_TV DESTINATION bin COMPONENT main)
