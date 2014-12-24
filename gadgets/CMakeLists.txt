if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETS__)
  link_directories(${Boost_LIBRARY_DIRS})
endif (WIN32)

# These dependencies should have been found already to even get in here, but we play safe
find_package(ACE)
find_package(FFTW3 COMPONENTS single double)
find_package(Ismrmrd)

include_directories(
  ${CMAKE_BINARY_DIR}/apps/gadgetron
  ${ACE_INCLUDE_DIR} 
  ${Boost_INCLUDE_DIR}
  ${FFTW3_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/apps/gadgetron 
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
  ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
  )

add_subdirectory(mri_core)
add_subdirectory(interventional_mri)
add_subdirectory(util)

if (ARMADILLO_FOUND)
    add_subdirectory(epi)
elseif (ARMADILLO_FOUND)
    MESSAGE("Armadillo not found, NOT compiling EPI Gadgets")
endif (ARMADILLO_FOUND)

if (MKL_FOUND OR ARMADILLO_FOUND)
    add_subdirectory(gtPlus)
endif (MKL_FOUND OR ARMADILLO_FOUND)

if (CUDA_FOUND)
  message("Cuda found, compiling gpu accelerated gadgets")
  add_subdirectory(pmri)
  add_subdirectory(grappa)
  add_subdirectory(radial)
  add_subdirectory(spiral)
else (CUDA_FOUND)
  message("Cuda NOT found, NOT compiling gpu accelerated gadgets")
endif(CUDA_FOUND)


if (MATLAB_FOUND)
  message("Compilng MATLAB gadgets")
  add_subdirectory(matlab)
endif(MATLAB_FOUND)

# nest the find_package checks to provide more useful error messages
find_package(Boost COMPONENTS python system thread REQUIRED)

if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    MESSAGE("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
    MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
    MESSAGE("NUMPY_INCLUDE_DIRS: ${NUMPY_INCLUDE_DIRS}")
    MESSAGE("Compiling Python Gadgets")
    add_subdirectory(python)
else (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    if(NOT PYTHONLIBS_FOUND)
        MESSAGE("Python Libraries/Headers NOT found, NOT compiling Python Gadgets")
    endif(NOT PYTHONLIBS_FOUND)
    if(NOT NUMPY_FOUND)
        MESSAGE("NumPy NOT found, NOT compiling Python Gadgets")
    endif(NOT NUMPY_FOUND)
    if(NOT Boost_PYTHON_FOUND)
        MESSAGE("Boost Python NOT found, NOT compiling Python Gadgets")
    endif(NOT Boost_PYTHON_FOUND)
endif (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)

find_package(DCMTK)
if(DCMTK_FOUND)
  message("Compiling DICOM gadget")
  add_subdirectory(dicom)
else(DCMTK_FOUND)
  message("DCMTK NOT found, not compiling DICOM gadget")
endif(DCMTK_FOUND)

add_subdirectory(cartesian)

if(ARMADILLO_FOUND)
    add_subdirectory(moco)
endif(ARMADILLO_FOUND)
