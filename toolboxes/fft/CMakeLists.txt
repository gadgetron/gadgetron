add_subdirectory(cpu)

if(CUDA_FOUND)
  add_subdirectory(gpu)
else()
  message("CUDA not found, not compiling gpu fft toolbox")
endif()