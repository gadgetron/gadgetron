if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    link_libraries(${MKL_LIBRARIES})
endif (MKL_FOUND)

if (ARMADILLO_FOUND)
  add_subdirectory(cpu)
endif (ARMADILLO_FOUND)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif(CUDA_FOUND)
