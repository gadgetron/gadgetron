#if(ACE_FOUND)
#  add_subdirectory(utilities)
#endif()

if(ISMRMRD_FOUND AND HDF5_FOUND)
  add_subdirectory(gadgetron_ismrmrd_client)
endif()