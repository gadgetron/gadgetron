include_directories(${CMAKE_SOURCE_DIR}/operators)
include_directories(${CMAKE_SOURCE_DIR}/solvers)

add_executable(radial_sense_sbc main_sbc.cpp)
add_executable(radial_sense_fista main_fista.cpp)

MESSAGE("CUDA LIBRARIES:  ${CUDA_LIBRARIES}")

target_link_libraries(radial_sense_sbc tomography_operators gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

target_link_libraries(radial_sense_fista tomography_operators gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})

