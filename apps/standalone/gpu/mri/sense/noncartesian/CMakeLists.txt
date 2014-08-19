include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  )

add_executable(sense_cg_generic_2d generic_cg.cpp)
target_link_libraries(sense_cg_generic_2d gpuoperators cpucore gpucore gpuparallelmri gpunfft hostutils gpusolvers ${CUDA_LIBRARIES})
install(TARGETS sense_cg_generic_2d DESTINATION bin)

add_subdirectory(radial)
