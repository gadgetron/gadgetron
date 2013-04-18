include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  )

add_subdirectory(2d)
add_subdirectory(3d)
