add_executable(register_CK_3d_cpu register_CK_3d.cpp)

target_link_libraries(register_CK_3d_cpu
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  armadillo
  )

add_executable(register_3d_volumes register_3d_volumes.cpp)

target_link_libraries(register_3d_volumes
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  armadillo
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${ARMADILLO_LIBRARIES} 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_mri_core
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpufft 
  gadgetron_toolbox_cpuklt 
  gadgetron_toolbox_cpuoperator 
  gadgetron_toolbox_cmr 
  gadgetron_toolbox_log
  )

  add_executable(applyDeformation applyDeformation.cpp)

target_link_libraries(applyDeformation
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_hostutils 
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_gpureg
  gadgetron_toolbox_gpuparallelmri
  armadillo
  ${ISMRMRD_LIBRARIES} 
  ${FFTW3_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${ARMADILLO_LIBRARIES} 
  gadgetron_toolbox_cpucore 
  gadgetron_toolbox_cpucore_math
  gadgetron_toolbox_mri_core
  gadgetron_toolbox_cpureg 
  gadgetron_toolbox_cpufft 
  gadgetron_toolbox_cpuklt 
  gadgetron_toolbox_cpuoperator 
  gadgetron_toolbox_cmr 
  gadgetron_toolbox_log
  )

install(TARGETS register_CK_3d_cpu DESTINATION bin COMPONENT main)
install(TARGETS register_3d_volumes DESTINATION bin COMPONENT main)

