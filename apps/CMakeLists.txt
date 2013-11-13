
IF (ACE_FOUND AND XSD_FOUND)
  add_subdirectory(gadgetron)
  IF (ISMRMRD_FOUND)
    add_subdirectory(clients)
  ENDIF (ISMRMRD_FOUND)	 
ENDIF (ACE_FOUND AND XSD_FOUND)
add_subdirectory(standalone)
