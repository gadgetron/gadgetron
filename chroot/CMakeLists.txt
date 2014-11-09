
configure_file("start-gadgetron.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/start-gadgetron.sh @ONLY)
configure_file("enter-chroot-env.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/enter-chroot-env.sh @ONLY)

configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_app.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web_app.cfg @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web.conf.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web.conf @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/apps/gadgetron/webapp/gadgetron_web_ld.conf.in" ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_web_ld.conf @ONLY)

configure_file("run-webapp.sh.in" ${CMAKE_CURRENT_BINARY_DIR}/run-webapp.sh @ONLY)


add_custom_target(chroot)

add_custom_command(
    TARGET chroot
    POST_BUILD
    VERBATIM
    COMMENT "Creating chroot tar file"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND sudo bash create_chroot.sh ${CMAKE_INSTALL_PREFIX}
                                       ${CMAKE_BINARY_DIR}
                                       ${CMAKE_SOURCE_DIR}
                                       ${GADGETRON_GIT_SHA1_HASH}
                                       ${ISMRMRD_LIB_DIR}
)

message("1. cmake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("2. cmake binary dir: " ${CMAKE_BINARY_DIR})
message("3. cmake source dir: " ${CMAKE_SOURCE_DIR})
message("ISMRMRD_HOME: " $ENV{ISMRMRD_HOME})
