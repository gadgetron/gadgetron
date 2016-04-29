IF (ACE_FOUND )
  add_subdirectory(gadgetron)
ENDIF (ACE_FOUND )

add_subdirectory(clients)
add_subdirectory(standalone)

IF (MATLAB_FOUND )
  add_subdirectory(matlab)
ENDIF (MATLAB_FOUND )
