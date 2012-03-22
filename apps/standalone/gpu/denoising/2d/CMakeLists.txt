if (WIN32)
ADD_DEFINITIONS(-D_USE_MATH_DEFINES)
endif (WIN32)

include_directories( ${CUDA_INCLUDE_DIRS}
                     ${Boost_INCLUDE_DIR}
                     ${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
                     ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
                     ${CMAKE_SOURCE_DIR}/toolboxes/solvers
                     ${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

add_executable(denoise_TV denoise_TV.cpp)

target_link_libraries(denoise_TV gpucore hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS denoise_TV DESTINATION bin)
