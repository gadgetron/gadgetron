if (CUDA_FOUND)
  add_subdirectory(gpu)
endif(CUDA_FOUND)

install (FILES 
  SenseJob.h
  DESTINATION include)
