if (CUDA_FOUND)
  add_subdirectory(gpu)
endif()

add_subdirectory(cpu)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)