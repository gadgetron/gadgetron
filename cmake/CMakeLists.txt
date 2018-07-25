install(FILES 	
    FindACE.cmake
    FindArmadillo.cmake
    FindCUDA_advanced.cmake
    FindDCMTK.cmake
    FindFFTW3.cmake
    FindGadgetron.cmake
    FindGLEW.cmake
    FindGMatlab.cmake
    FindGperftools.cmake
    FindMKL.cmake
    FindNumPy.cmake
    FindOctave.cmake
    FindPLplot.cmake
    FindPython3.cmake
    FindZFP.cmake
    DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)

if (WIN32)
    install(FILES InstallWinGadgetron.bat DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)
endif ()

#if (UNIX)
#    if (NOT APPLE)
#        install(FILES ./debian/postinst DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT main)
#        install(FILES ./debian/postinst_web DESTINATION ${GADGETRON_INSTALL_CMAKE_PATH} COMPONENT web)
#    endif ()
#endif ()
