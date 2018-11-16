if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()

add_subdirectory(cpu)
