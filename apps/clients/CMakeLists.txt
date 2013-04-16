if(WIN32)
  add_definitions(-DTIXML_USE_STL)
endif(WIN32)
 
if(ISMRMRD_FOUND AND HDF5_FOUND)
  add_subdirectory(mriclient)
endif(ISMRMRD_FOUND AND HDF5_FOUND)
