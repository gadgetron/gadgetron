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
  message("Cuda found, compiling gpu accelerated gadgets")
  add_subdirectory(radial)
  add_subdirectory(spiral)
  add_subdirectory(sense)
  add_subdirectory(grappa)
else (CUDA_FOUND)
  message("Cuda NOT found, NOT compiling gpu accelerated gadgets")
endif(CUDA_FOUND)

#find_package(Octave)
#if (OCTAVE_FOUND)
#  add_subdirectory(octave)
#endif(OCTAVE_FOUND)

if (MATLAB_FOUND)
  message("Compilng MATLAB gadgets")
  add_subdirectory(matlab)
endif(MATLAB_FOUND)

# nest the find_package checks to provide more useful error messages
find_package(Boost COMPONENTS python)

if (Boost_PYTHON_FOUND)
    if (PYTHONLIBS_FOUND)
        MESSAGE("PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
        MESSAGE("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
        if (NUMPY_FOUND)
            MESSAGE("NUMPY_INCLUDE_DIRS: ${NUMPY_INCLUDE_DIRS}")
            if (NUMPY_VERSION VERSION_LESS "1.7")
                message("NumPy version < 1.7. NOT compiling Python Gadgets")
            else (${NUMPY_VERSION} VERSION_LESS "1.7")
                MESSAGE("Compiling Python Gadgets")
                add_subdirectory(python)
            endif (NUMPY_VERSION VERSION_LESS "1.7")
        else(NUMPY_FOUND)
            MESSAGE("NumPy NOT found, NOT compiling Python Gadgets")
        endif(NUMPY_FOUND)
    else (PYTHONLIBS_FOUND)
        MESSAGE("Python Libraries/Headers NOT found, NOT compiling Python Gadgets")
    endif (PYTHONLIBS_FOUND)
else (Boost_PYTHON_FOUND)
    MESSAGE("Boost Python NOT found, NOT compiling Python Gadgets")
endif (Boost_PYTHON_FOUND)

find_package(DCMTK)
if(DCMTK_FOUND)
  message("Compiling DICOM gadget")
  add_subdirectory(dicom)
else(DCMTK_FOUND)
  message("DCMTK NOT found, not compiling DICOM gadget")
endif(DCMTK_FOUND)

add_subdirectory(cartesian)
add_subdirectory(moco)
