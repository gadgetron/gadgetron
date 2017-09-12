include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
  ${ACE_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  )

if(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.819" )
add_subdirectory(2d)
add_subdirectory(3d)
endif()
