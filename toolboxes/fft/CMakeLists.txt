add_subdirectory(cpu)

add_compile_options(-Werror)
if(CUDA_FOUND)
  add_subdirectory(gpu)
else()
  message("CUDA not found, not compiling gpu fft toolbox")
endif()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)