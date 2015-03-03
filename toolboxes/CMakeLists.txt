if (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
endif (WIN32)

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    link_libraries(${MKL_LIBRARIES})
endif (MKL_FOUND)

add_subdirectory(log)
add_subdirectory(core)
add_subdirectory(mri_core)

add_subdirectory(operators)
add_subdirectory(solvers)

add_subdirectory(fft)
add_subdirectory(nfft)

add_subdirectory(dwt)

add_subdirectory(mri)
add_subdirectory(ct)

add_subdirectory(registration)

IF (ACE_FOUND)
  add_subdirectory(gadgettools)
  add_subdirectory(cloudbus)
ENDIF()

IF (FFTW3_FOUND AND ISMRMRD_FOUND)
  add_subdirectory(gtplus)
ENDIF()

find_package(BLAS)
find_package(LAPACK)

# Should we compile the python toolbox
find_package(Boost COMPONENTS python system thread REQUIRED)
if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    add_subdirectory(python)
endif (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
