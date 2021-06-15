set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "")

find_package(CUDA)

# Check for GPUs present and their compute capability
# based on http://stackoverflow.com/questions/2285185/easiest-way-to-test-for-existence-of-cuda-capable-gpu-from-cmake/2297877#2297877 (Christopher Bruns)
if(CUDA_FOUND)
    if (${CUDA_VERSION_MAJOR} VERSION_GREATER "7")
        set(CUDA_NVCC_FLAGS6 "-gencode arch=compute_60,code=sm_60")   
        set(CUDA_NVCC_FLAGS61 "-gencode arch=compute_61,code=sm_61")   
    endif()
    if (${CUDA_VERSION_MAJOR} VERSION_GREATER "8")
        set(CUDA_NVCC_FLAGS7 "-gencode arch=compute_70,code=sm_70")
    endif()
        if (${CUDA_VERSION_MAJOR} VERSION_GREATER "9")
                set(CUDA_NVCC_FLAGS8 "-gencode arch=compute_80,code=sm_80")

    endif()

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
    set(CUDA_COMPUTE_CAPABILITY ${RUN_OUTPUT_VAR} CACHE STRING "Compute capability of CUDA-capable GPU present. Seperate multiple by ;. For all known, use ALL")
  else()
    set(CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present")
    set(CUDA_COMPUTE_CAPABILITY ALL CACHE STRING "Compute capability of CUDA-capable GPU present. Seperate multiple by ;. For all known, use ALL")
  endif()

find_cuda_helper_libs(cusparse)
set(CUDA_CUSPARSE_LIBRARIES ${CUDA_cusparse_LIBRARY})
if( "${CUDA_COMPUTE_CAPABILITY}" MATCHES ALL)
  if (${CUDA_VERSION_MAJOR} VERSION_GREATER "6")
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS4} ${CUDA_NVCC_FLAGS5} ${CUDA_NVCC_FLAGS52} ${CUDA_NVCC_FLAGS6} ${CUDA_NVCC_FLAGS61} ${CUDA_NVCC_FLAGS7} ${CUDA_NVCC_FLAGS8})
  else()
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_NVCC_FLAGS4})
  endif()
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

