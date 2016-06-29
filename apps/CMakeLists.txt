if (ACE_FOUND )
  add_subdirectory(gadgetron)
endif ()

add_subdirectory(clients)
add_subdirectory(standalone)

if (MATLAB_FOUND )
  add_subdirectory(matlab)
endif ()
