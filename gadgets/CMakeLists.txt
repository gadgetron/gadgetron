add_subdirectory(core)

if (CUDA_FOUND)
add_subdirectory(grappa)
add_subdirectory(cgsense)
if (NOT WIN32)
	add_subdirectory(spiral)
endif (NOT WIN32)
endif(CUDA_FOUND)

if (BLAS_FOUND AND LAPACK_FOUND)
	MESSAGE("Compiling optimized (BLAS/LAPACK) MRI Gadgets")
	add_subdirectory(mri_core)
endif (BLAS_FOUND AND LAPACK_FOUND)

if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND)
   MESSAGE("Compiling Python Gadgets")
   add_subdirectory(python)
endif(Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND)
