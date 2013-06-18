if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
  link_directories(${Boost_LIBRARY_DIRS})
endif (WIN32)

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

add_subdirectory(mri_core)

if (CUDA_FOUND)
  add_subdirectory(radial)
  add_subdirectory(spiral)
  add_subdirectory(sense)
  add_subdirectory(grappa)
endif(CUDA_FOUND)

#find_package(Octave)
#if (OCTAVE_FOUND)
#  add_subdirectory(octave)
#endif(OCTAVE_FOUND)

if (MATLAB_FOUND)
  message("MATLAB FOUND: ${MATLAB_INCLUDE_DIR}")
  add_subdirectory(matlab)
else(MATLAB_FOUND)
  message("MATLAB NOT FOUND: ${MATLAB_INCLUDE_DIR}")
endif(MATLAB_FOUND)

find_package(Boost COMPONENTS python)
if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND PYTHON_NUMPY_FOUND)
  MESSAGE("Compiling Python Gadgets")
  MESSAGE("PYTHON_NUMPY_INCLUDE_DIR: ${PYTHON_NUMPY_INCLUDE_DIR}")
  MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
  MESSAGE("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
  add_subdirectory(python)
endif(Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND PYTHON_NUMPY_FOUND)
