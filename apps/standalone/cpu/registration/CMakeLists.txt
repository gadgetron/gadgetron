include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  )

if(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.79" )
add_subdirectory(2d)
add_subdirectory(3d)
endif(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.79" )
