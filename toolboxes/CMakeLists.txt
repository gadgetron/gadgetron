if (WIN32)
    add_definitions(-DTIXML_USE_STL)
endif ()

if (MKL_FOUND)
    include_directories( ${MKL_INCLUDE_DIR} )
    link_directories( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    link_libraries(${MKL_LIBRARIES})
endif ()

add_subdirectory(log)

add_subdirectory(operators)
add_subdirectory(solvers)

if (FFTW3_FOUND)
  add_subdirectory(fft)
  add_subdirectory(core)
  if(ARMADILLO_FOUND)
    add_subdirectory(mri_core)
  endif()
  add_subdirectory(mri)
  if(ARMADILLO_FOUND)
    add_subdirectory(cmr)
  endif()
  if (NOT(WIN32))
    add_subdirectory(ct)
  endif()

  if (ISMRMRD_FOUND)
    add_subdirectory(mri_image)
  endif()
endif()

if(ARMADILLO_FOUND)
  add_subdirectory(klt)
endif()
add_subdirectory(nfft)
add_subdirectory(dwt)
add_subdirectory(registration)
add_subdirectory(rest)
add_subdirectory(fatwater)
add_subdirectory(ffd)
add_subdirectory(image)
add_subdirectory(pattern_recognition)

if (ACE_FOUND)
  add_subdirectory(gadgettools)
  add_subdirectory(cloudbus)
endif()

add_subdirectory(image_io)

find_package(BLAS)
find_package(LAPACK)


if (MATLAB_FOUND)
  message("Compiling MATLAB toolbox")
  add_subdirectory(matlab)
endif()


# Should we compile the python toolbox
find_package(Boost COMPONENTS python system thread REQUIRED)
if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    add_subdirectory(python)
endif ()

if (ARMADILLO_FOUND AND Boost_FOUND AND PLPLOT_FOUND)
    add_subdirectory(plplot)
endif ()
