add_subdirectory(core)

if (CUDA_FOUND)
add_subdirectory(grappa)
add_subdirectory(cgsense)
endif(CUDA_FOUND)

if (NOT WIN32)
add_subdirectory(spiral)
endif (NOT WIN32)

if (PYTHONLIBS_FOUND AND NOT WIN32)
   add_subdirectory(python)
endif(PYTHONLIBS_FOUND AND NOT WIN32)
