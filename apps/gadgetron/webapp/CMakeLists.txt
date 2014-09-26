
configure_file("${CMAKE_SOURCE_DIR}/cmake/gadgetron_web_app.in" ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/cmake/gadgetron_web.conf.in" ${CMAKE_BINARY_DIR}/gadgetron_web.conf @ONLY)

if (WIN32)
    install(FILES gadgetron_web_app.py DESTINATION bin)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION config)
else (WIN32)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web.conf DESTINATION /etc/init)
    install(FILES gadgetron_web_app.py DESTINATION bin)
    install(FILES ${CMAKE_BINARY_DIR}/gadgetron_web_app.cfg DESTINATION config)
endif (WIN32)
