
configure_file("start-gadgetron.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron.sh @ONLY)
configure_file("enter-chroot-env.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/enter-chroot-env.sh @ONLY)

add_custom_target(chroot)

add_custom_command(
    TARGET chroot
    POST_BUILD
    VERBATIM
    COMMENT "Creating chroot tar file"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND sh create_chroot.sh ${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}
)

message("cmake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("cmake current binary dir: " ${CMAKE_CURRENT_BINARY_DIR})
message("cmake current source dir: " ${CMAKE_CURRENT_SOURCE_DIR})
message("cmake binary dir: " ${CMAKE_BINARY_DIR})
message("cmake source dir: " ${CMAKE_SOURCE_DIR})
message("concatenation: " ${CMAKE_CURRENT_SOURCE_DIR}/chroot-root${CMAKE_INSTALL_PREFIX})
