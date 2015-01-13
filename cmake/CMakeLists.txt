install(FILES 	
FindACE.cmake
FindFFTW3.cmake
FindGLEW.cmake
FindNumPy.cmake
FindGadgetron.cmake
FindIsmrmrd.cmake
FindGMatlab.cmake
FindDCMTK.cmake
FindMKL.cmake
DESTINATION cmake COMPONENT main)

if (WIN32)
    install(FILES InstallWinGadgetron.bat DESTINATION cmake COMPONENT main)
endif (WIN32)

#if (UNIX)
#    if (NOT APPLE)
#        install(FILES ./debian/postinst DESTINATION cmake COMPONENT main)
#        install(FILES ./debian/postinst_web DESTINATION cmake COMPONENT web)
#    endif (NOT APPLE)
#endif (UNIX)
