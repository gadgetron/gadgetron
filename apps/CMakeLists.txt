set(CMAKE_BUILD_RPATH ${GADGETRON_BUILD_RPATH})

add_subdirectory(gadgetron)
add_subdirectory(storage)

add_subdirectory(clients)
add_subdirectory(standalone)
add_subdirectory(matlab)

