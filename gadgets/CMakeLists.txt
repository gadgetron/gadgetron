add_subdirectory(core)

if (CUDA_FOUND)
add_subdirectory(grappa)
add_subdirectory(cgsense)
endif(CUDA_FOUND)

if (NOT WIN32)
add_subdirectory(spiral)
endif (NOT WIN32)

if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND)
   MESSAGE("Compiling Python Gadgets")
   add_subdirectory(python)
endif(Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND)
