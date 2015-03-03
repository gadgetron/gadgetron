IF(ACE_FOUND)
  add_subdirectory(utilities)
ENDIF()

IF(ISMRMRD_FOUND AND ACE_FOUND AND HDF5_FOUND)
  add_subdirectory(gadgetron_ismrmrd_client)
ENDIF()