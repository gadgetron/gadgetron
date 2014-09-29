install(FILES 	
FindACE.cmake
FindFFTW3.cmake
FindGLEW.cmake
FindNumPy.cmake
FindTinyXML.cmake
FindGadgetron.cmake
FindIsmrmrd.cmake
FindGMatlab.cmake
FindDCMTK.cmake
FindMKL.cmake
DESTINATION cmake)

if (WIN32)
    install(FILES InstallWinGadgetron.bat DESTINATION cmake)
endif (WIN32)
