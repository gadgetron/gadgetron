add_subdirectory(core)

if (CUDA_FOUND)
add_subdirectory(grappa)
add_subdirectory(cgsense)
add_subdirectory(spiral)
endif(CUDA_FOUND)

if (PYTHONLIBS_FOUND)
   add_subdirectory(python)
endif(PYTHONLIBS_FOUND)
