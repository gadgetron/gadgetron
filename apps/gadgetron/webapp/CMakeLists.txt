configure_file("gadgetron_web_app.in" ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg @ONLY)
configure_file("gadgetron_web.conf.in" ${CMAKE_BINARY_DIR}/gadgetron_web.conf @ONLY)
configure_file("gadgetron_web_ld.conf.in" ${CMAKE_BINARY_DIR}/gadgetron_web_ld.conf @ONLY)

if (WIN32)
    install(FILES gadgetron_web_app.py DESTINATION bin)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH})
else ()
    #install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web.conf DESTINATION /etc/init COMPONENT web)
    #install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_ld.conf DESTINATION /etc/ld.so.conf.d COMPONENT web)
    install(FILES gadgetron_web_app.py DESTINATION bin COMPONENT web)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION ${GADGETRON_INSTALL_CONFIG_PATH} COMPONENT web)
endif ()

message("Add gadgetron_web_info ...")
add_executable(gadgetron_web_info main.cpp)
install(TARGETS gadgetron_web_info DESTINATION bin COMPONENT web)

install(FILES
    ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg
    ${CMAKE_BINARY_DIR}/gadgetron_web.conf
    ${CMAKE_BINARY_DIR}/gadgetron_web_ld.conf
    ${CMAKE_SOURCE_DIR}/chroot/gadgetron_chroot.conf
    gadgetron_web_app.py
    DESTINATION ${GADGETRON_INSTALL_CHROOT_SCRIPTS_PATH}
    COMPONENT scripts
    )
