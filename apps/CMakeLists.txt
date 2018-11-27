add_subdirectory(gadgetron)

add_subdirectory(clients)
add_subdirectory(standalone)

if (MATLAB_FOUND )
  add_subdirectory(matlab)
endif ()
