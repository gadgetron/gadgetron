if (WIN32)
    ADD_DEFINITIONS(-DTIXML_USE_STL)
endif (WIN32)

if (MKL_FOUND)
    INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
    LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    link_libraries(${MKL_LIBRARIES})
endif (MKL_FOUND)

add_subdirectory(log)

add_subdirectory(operators)
add_subdirectory(solvers)

if (FFTW3_FOUND)
  add_subdirectory(fft)
  add_subdirectory(core)
  add_subdirectory(mri_core)
  add_subdirectory(mri)
  add_subdirectory(cmr)
  if (NOT(WIN32))
    add_subdirectory(ct)
  endif()

  if (ISMRMRD_FOUND)
    add_subdirectory(mri_image)
  endif(ISMRMRD_FOUND)
endif()

add_subdirectory(klt)
add_subdirectory(nfft)
add_subdirectory(dwt)
add_subdirectory(registration)
add_subdirectory(rest)
add_subdirectory(fatwater)
add_subdirectory(ffd)

IF (ACE_FOUND)
  add_subdirectory(gadgettools)
  add_subdirectory(cloudbus)
ENDIF()

add_subdirectory(image_io)

find_package(BLAS)
find_package(LAPACK)


if (MATLAB_FOUND)
  message("Compiling MATLAB toolbox")
  add_subdirectory(matlab)
endif(MATLAB_FOUND)


# Should we compile the python toolbox
find_package(Boost COMPONENTS python system thread REQUIRED)
if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    add_subdirectory(python)
endif (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)

if (ARMADILLO_FOUND AND Boost_FOUND AND PLPLOT_FOUND)
    add_subdirectory(plplot)
endif (ARMADILLO_FOUND AND Boost_FOUND AND PLPLOT_FOUND)
