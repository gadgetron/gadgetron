install(FILES 	
FindACE.cmake
FindFFTW3.cmake
FindGLEW.cmake
FindNumPy.cmake
FindGadgetron.cmake
FindGMatlab.cmake
FindDCMTK.cmake
FindMKL.cmake
FindPLplot.cmake
DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)

if (WIN32)
    install(FILES InstallWinGadgetron.bat DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)
endif (WIN32)

#if (UNIX)
#    if (NOT APPLE)
#        install(FILES ./debian/postinst DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)
#        install(FILES ./debian/postinst_web DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT web)
#    endif (NOT APPLE)
#endif (UNIX)
