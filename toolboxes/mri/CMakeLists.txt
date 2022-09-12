add_subdirectory(pmri)
add_subdirectory(spiral)

add_subdirectory(epi)

add_subdirectory(sdc)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)