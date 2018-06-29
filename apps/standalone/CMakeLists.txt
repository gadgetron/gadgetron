if (MKL_FOUND)
    include_directories( ${MKL_INCLUDE_DIR} )
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    link_libraries(${MKL_LIBRARIES})
endif ()
include_directories(<
    ${CMAKE_BINARY_DIR}/toolboxes/core
)

if (ARMADILLO_FOUND)
  add_subdirectory(cpu)
endif ()

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif()
