if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

add_executable(denoise_TV denoise_TV.cpp)

target_link_libraries(denoise_TV gadgetron_toolbox_gpucore gadgetron_toolbox_hostutils gadgetron_toolbox_gpusolvers gadgetron_toolbox_gpuoperators ${CUDA_LIBRARIES})

install(TARGETS denoise_TV DESTINATION bin)
