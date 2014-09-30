if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

add_executable(denoise_TV denoise_TV.cpp)

target_link_libraries(denoise_TV gpucore hostutils gpusolvers gpuoperators ${CUDA_LIBRARIES})

install(TARGETS denoise_TV DESTINATION bin COMPONENT main)
