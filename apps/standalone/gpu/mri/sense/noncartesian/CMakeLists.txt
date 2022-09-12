include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
  )

add_executable(sense_cg_generic_2d generic_cg.cpp)

target_link_libraries(sense_cg_generic_2d gadgetron_toolbox_gpu gadgetron_toolbox_cpucore gadgetron_toolbox_hostutils ${CUDA_LIBRARIES})

install(TARGETS sense_cg_generic_2d DESTINATION bin COMPONENT main)

add_subdirectory(radial)
