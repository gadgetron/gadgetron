include_directories(   
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${ISMRMRD_XSD_INCLUDE_DIR}
  ${XSD_INCLUDE_DIR}
  ${XERCESC_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
  ${CMAKE_SOURCE_DIR}/apps/gadgetron 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  )

add_subdirectory(core)

if (CUDA_FOUND)
  add_subdirectory(grappa)
  add_subdirectory(cgsense)
  add_subdirectory(spiral)
endif(CUDA_FOUND)

#find_package(Octave)
#if (OCTAVE_FOUND)
#  add_subdirectory(octave)
#endif(OCTAVE_FOUND)

find_package(GMatlab)
if (MATLAB_FOUND)
  message("MATLAB FOUND: ${MATLAB_INCLUDE_DIR}")
  add_subdirectory(matlab)
else(MATLAB_FOUND)
  message("MATLAB NOT FOUND: ${MATLAB_INCLUDE_DIR}")
endif(MATLAB_FOUND)

if (ARMADILLO_FOUND)
  MESSAGE("Compiling optimized (Armadillo) MRI Gadgets")
  add_subdirectory(mri_core)
endif (ARMADILLO_FOUND)

find_package(Boost COMPONENTS python)
if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND PYTHON_NUMPY_FOUND)
  MESSAGE("Compiling Python Gadgets")
  MESSAGE("PYTHON_NUMPY_INCLUDE_DIR: ${PYTHON_NUMPY_INCLUDE_DIR}")
  MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
  MESSAGE("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
  add_subdirectory(python)
endif(Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND PYTHON_NUMPY_FOUND)
