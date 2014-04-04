IF (ACE_FOUND AND XSD_FOUND AND XERCESC_FOUND )
  add_subdirectory(gadgetron)
ENDIF (ACE_FOUND AND XSD_FOUND AND XERCESC_FOUND)

add_subdirectory(clients)
add_subdirectory(standalone)
