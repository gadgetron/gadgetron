add_subdirectory(cpucore)
add_subdirectory(hostutils)

IF (ACE_FOUND)
        add_subdirectory(gadgettools)
ENDIF (ACE_FOUND)

if (BLAS_FOUND)
 	add_subdirectory(linalg)
endif(BLAS_FOUND)

if (CUDA_FOUND)
   add_subdirectory(gpucore)
   add_subdirectory(gpunfft)
   add_subdirectory(gpupmri)
   add_subdirectory(solvers)
endif(CUDA_FOUND)

FIND_PACKAGE(HDF5 1.8 COMPONENTS C CXX)

if (HDF5_FOUND)
   MESSAGE("-- HDF5 Support found, Compiling HDF5 Support") 
   add_subdirectory(hdf5utils)
   add_subdirectory(hdf5mri)
endif(HDF5_FOUND)

