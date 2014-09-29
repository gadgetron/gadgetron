
configure_file("gadgetron_web_app.in" ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg @ONLY)
configure_file("gadgetron_web.conf.in" ${CMAKE_BINARY_DIR}/gadgetron_web.conf @ONLY)
configure_file("gadgetron_web_ld.conf.in" ${CMAKE_BINARY_DIR}/gadgetron_web_ld.conf @ONLY)

if (WIN32)
    install(FILES gadgetron_web_app.py DESTINATION bin)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION config)
else (WIN32)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web.conf DESTINATION /etc/init)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_ld.conf DESTINATION /etc/ld.so.conf.d)
    install(FILES gadgetron_web_app.py DESTINATION bin)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION config)
endif (WIN32)
