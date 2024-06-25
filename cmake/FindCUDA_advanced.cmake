set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "")

find_package(CUDA 12.3)

# Check for GPUs present and their compute capability
# based on http://stackoverflow.com/questions/2285185/easiest-way-to-test-for-existence-of-cuda-capable-gpu-from-cmake/2297877#2297877 (Christopher Bruns)
if(CUDA_FOUND)

  # Enumerate the compute capabilities we will be building for if not targeting build system GPU
  set(CUDA_NVCC_FLAGS7  "-gencode arch=compute_70,code=sm_70")
  set(CUDA_NVCC_FLAGS75 "-gencode arch=compute_75,code=sm_75")
  set(CUDA_NVCC_FLAGS8  "-gencode arch=compute_80,code=sm_80")
  set(CUDA_NVCC_FLAGS86 "-gencode arch=compute_86,code=sm_86")
  set(CUDA_NVCC_FLAGS87 "-gencode arch=compute_87,code=sm_87")
  set(CUDA_NVCC_FLAGS90 "-gencode arch=compute_90,code=sm_90")

  cuda_find_helper_file(cuda_compute_capability c)
  try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
    ${CMAKE_BINARY_DIR} 
    ${CUDA_cuda_compute_capability}
    CMAKE_FLAGS 
    -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
    -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
    COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
    RUN_OUTPUT_VARIABLE RUN_OUTPUT_VAR)
  # COMPILE_RESULT_VAR is TRUE when compile succeeds
  # RUN_RESULT_VAR is zero when a GPU is found
  if(COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR)
    set(CUDA_HAVE_GPU TRUE CACHE BOOL "Whether CUDA-capable GPU is present")
    set(CUDA_COMPUTE_CAPABILITY ${RUN_OUTPUT_VAR} CACHE STRING "Compute capability of CUDA-capable GPU present. Separate multiple by ;. For all known, use ALL")
  else()
    set(CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present")
    set(CUDA_COMPUTE_CAPABILITY ALL CACHE STRING "Compute capability of CUDA-capable GPU present. Separate multiple by ;. For all known, use ALL")
  endif()

find_cuda_helper_libs(cusparse)
set(CUDA_CUSPARSE_LIBRARIES ${CUDA_cusparse_LIBRARY})
if( "${CUDA_COMPUTE_CAPABILITY}" MATCHES ALL)
  set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS} 
    ${CUDA_NVCC_FLAGS7}
    ${CUDA_NVCC_FLAGS75}
    ${CUDA_NVCC_FLAGS8}
    ${CUDA_NVCC_FLAGS86}
    ${CUDA_NVCC_FLAGS87}
    ${CUDA_NVCC_FLAGS90})
else()
  set(CUDA_MOSTUSED_ARCH "")
  foreach(code ${CUDA_COMPUTE_CAPABILITY})
    if (NOT CUDA_MOSTUSED_ARCH)
      set (CUDA_MOSTUSED_ARCH "-arch=sm_${code}")
      set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_${code}")
    endif()
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_${code},code=sm_${code} ")
  endforeach()
endif()

endif()
