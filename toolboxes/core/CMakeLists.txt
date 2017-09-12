include_directories(
  ${Boost_INCLUDE_DIR}
  )

configure_file(core_defines.h.in ${CMAKE_CURRENT_SOURCE_DIR}/core_defines.h)

install(FILES 
  core_defines.h
  NDArray.h
  complext.h
  vector_td.h
  vector_td_operators.h
  vector_td_utilities.h
  vector_td_io.h
  real_utilities.h
  GadgetronException.h
  GadgetronTimer.h
  Gadgetron_enable_types.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(cpu)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()
