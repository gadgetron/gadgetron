add_subdirectory(pmri)

add_subdirectory(epi)
if (CUDA_FOUND)
    add_subdirectory(hyper)
endif ()

