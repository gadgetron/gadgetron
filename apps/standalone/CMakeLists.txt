include_directories(
    ${CMAKE_BINARY_DIR}/toolboxes/core
)


add_subdirectory(cpu)
if (CUDA_FOUND)
  add_subdirectory(gpu)
endif()
