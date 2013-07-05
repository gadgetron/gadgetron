if (ARMADILLO_FOUND)
  add_subdirectory(cpu)
endif (ARMADILLO_FOUND)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif(CUDA_FOUND)
