IF (WIN32)
ADD_DEFINITIONS(-DTIXML_USE_STL)
ENDIF (WIN32)
 
add_subdirectory(mriclient)
