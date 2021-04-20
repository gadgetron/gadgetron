include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  )

add_executable(sense_cg_generic_2d generic_cg.cpp)

target_link_libraries(sense_cg_generic_2d gadgetron_toolbox_gpuoperators gadgetron_toolbox_cpucore gadgetron_toolbox_gpucore gadgetron_toolbox_gpuparallelmri gadgetron_toolbox_gpunfft gadgetron_toolbox_hostutils gadgetron_toolbox_gpusolvers ${CUDA_LIBRARIES} ${PYTHON_LIBRARIES})

install(TARGETS sense_cg_generic_2d DESTINATION bin COMPONENT main)

add_subdirectory(radial)
