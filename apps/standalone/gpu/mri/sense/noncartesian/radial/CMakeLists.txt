include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/pmri/gpu
)

add_subdirectory(2d_golden_ratio)
add_subdirectory(2d_golden_ratio_kt)

#if (QT4_FOUND AND GLEW_FOUND AND GLUT_FOUND AND OPENGL_FOUND)
#  add_subdirectory(2d_golden_ratio_gui)
#endif ()
