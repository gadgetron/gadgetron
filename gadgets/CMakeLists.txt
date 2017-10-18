if (WIN32)
  add_definitions(-D__BUILD_GADGETS__)
  link_directories(${Boost_LIBRARY_DIRS})
endif ()

# These dependencies should have been found already to even get in here, but we play safe
find_package(ACE)
find_package(FFTW3 COMPONENTS single double)

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
add_subdirectory(fatwater)
add_subdirectory(interventional_mri)
if (ARMADILLO_FOUND)
  add_subdirectory(cmr)
endif()

find_program(BARTCMD bart)
if (BARTCMD)
	add_subdirectory(bart)
endif()

if (ARMADILLO_FOUND)
    add_subdirectory(epi)
elseif (ARMADILLO_FOUND)
    message("Armadillo not found, NOT compiling EPI Gadgets")
endif ()

if (CUDA_FOUND)
  message("Cuda found, compiling gpu accelerated gadgets")
  add_subdirectory(pmri)
  add_subdirectory(radial)  
  add_subdirectory(hyper)
  add_subdirectory(gpu)
  add_subdirectory(mri_noncartesian)
else ()
  message("Cuda NOT found, NOT compiling gpu accelerated gadgets")
endif()

add_subdirectory(spiral)
add_subdirectory(grappa)
add_subdirectory(distributed)

if (MATLAB_FOUND)
  message("Compiling MATLAB gadgets")
  add_subdirectory(matlab)
endif()

# nest the find_package checks to provide more useful error messages
find_package(Boost COMPONENTS python system thread REQUIRED)

if (Boost_PYTHON_FOUND AND PYTHONLIBS_FOUND AND NUMPY_FOUND)
    message("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
    message("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
    message("NUMPY_INCLUDE_DIRS: ${NUMPY_INCLUDE_DIRS}")
    message("Compiling Python Gadgets")
    add_subdirectory(python)   
    add_definitions(-DCOMPILING_WITH_PYTHON_SUPPORT)
else ()
    if(NOT PYTHONLIBS_FOUND)
        message("Python Libraries/Headers NOT found, NOT compiling Python Gadgets")
    endif()
    if(NOT NUMPY_FOUND)
        message("NumPy NOT found, NOT compiling Python Gadgets")
    endif()
    if(NOT Boost_PYTHON_FOUND)
        message("Boost Python NOT found, NOT compiling Python Gadgets")
    endif()
endif ()

add_subdirectory(util)

find_package(DCMTK)
if(DCMTK_FOUND)
  message("Compiling DICOM gadget")
  add_subdirectory(dicom)
else()
  message("DCMTK NOT found, not compiling DICOM gadget")
endif()

add_subdirectory(cartesian)

if(ARMADILLO_FOUND)
    add_subdirectory(moco)
endif()

if(PLPLOT_FOUND)
    if (MKL_FOUND OR ARMADILLO_FOUND)
        add_subdirectory(plplot)
    endif ()
endif()
