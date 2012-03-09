IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)
 
add_subdirectory(mriclient)

IF (QT4_FOUND)
add_subdirectory(imageclient)
ENDIF (QT4_FOUND)
