# - Try to find GLEW
# Once done this will define
#  
#  GLEW_FOUND        - system has GLEW
#  GLEW_INCLUDE_DIR  - the GLEW include directory
#  GLEW_LIBRARY_DIR  - where the libraries are
#  GLEW_LIBRARY      - Link these to use GLEW
#   

if (GLEW_INCLUDE_DIR)
  # Already in cache, be silent
  set(GLEW_FIND_QUIETLY TRUE)
endif ()

if( WIN32 )
   if( MSVC80 )
       set( COMPILER_PATH "C:/Program\ Files/Microsoft\ Visual\ Studio\ 8/VC" )
   endif()
   if( MSVC71 )
       set( COMPILER_PATH "C:/Program\ Files/Microsoft\ Visual\ Studio\ .NET\ 2003/Vc7" )
   endif()
   find_path( GLEW_INCLUDE_DIR gl/glew.h gl/wglew.h
              PATHS c:/glew/include ${COMPILER_PATH}/PlatformSDK/Include )
   set( GLEW_NAMES glew32 )
   find_library( GLEW_LIBRARY
                 NAMES ${GLEW_NAMES}
                 PATHS c:/glew/lib ${COMPILER_PATH}/PlatformSDK/Lib )
else()
   find_path( GLEW_INCLUDE_DIR glew.h wglew.h
              PATHS /usr/local/include /usr/include
              PATH_SUFFIXES gl/ GL/ )
   set( GLEW_NAMES glew GLEW )
   find_library( GLEW_LIBRARY
                 NAMES ${GLEW_NAMES}
                 PATHS /usr/lib /usr/local/lib )
endif()

get_filename_component( GLEW_LIBRARY_DIR ${GLEW_LIBRARY} PATH )

if (GLEW_INCLUDE_DIR AND GLEW_LIBRARY)
   set(GLEW_FOUND TRUE)
    set( GLEW_LIBRARY_DIR ${GLEW_LIBRARY} )
    message("GLEW FOUND")
else ()
   set( GLEW_FOUND FALSE )
   set( GLEW_LIBRARY_DIR )
    message("GLEW NOT FOUND")
endif ()

mark_as_advanced(
  GLEW_LIBRARY
  GLEW_INCLUDE_DIR
)
