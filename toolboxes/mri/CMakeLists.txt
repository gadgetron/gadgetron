add_subdirectory(pmri)

add_subdirectory(epi)
if (CUDA_FOUND)
    add_subdirectory(hyper)
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)