# - Find CULA
# Find the native CULA includes and library
#
#   CULA_FOUND       - True if CULA found.
#   CULA_INCLUDE_DIR - where to find cula.h, etc.
#   CULA_LIBRARIES   - List of libraries when using TinyXML.
#

IF( CULA_INCLUDE_DIR )
    # Already in cache, be silent
    SET( CULA_FIND_QUIETLY TRUE )
ENDIF( CULA_INCLUDE_DIR )

FIND_PATH( CULA_INCLUDE_DIR "cula.h"
           PATH_SUFFIXES "cula/include" )

MESSAGE("CULA_INCLUDE_DIR = ${CULA_INCLUDE_DIR}")


FIND_LIBRARY( CULA_LIBRARY
              NAMES "cula"
              PATH_SUFFIXES "cula/lib64" )

FIND_LIBRARY( CULA_LAPACK_LIBRARY
              NAMES "cula_lapack"
              PATH_SUFFIXES "cula/lib64" )

FIND_LIBRARY( CULA_CORE_LIBRARY
              NAMES "cula_core"
              PATH_SUFFIXES "cula/lib64" )

#This is version 12 of CULA
if (CULA_LIBRARY)
  list(APPEND CULA_LIBRARIES ${CULA_LIBRARY})
endif (CULA_LIBRARY)

#This is version 13 of CULA
if (CULA_LAPACK_LIBRARY)
  list(APPEND CULA_LIBRARIES ${CULA_LAPACK_LIBRARY})
endif (CULA_LAPACK_LIBRARY)

#This is version 13 of CULA
if (CULA_CORE_LIBRARY)
  list(APPEND CULA_LIBRARIES ${CULA_CORE_LIBRARY})
endif (CULA_CORE_LIBRARY)

MESSAGE("CULA_LIBRARIES = ${CULA_LIBRARIES}")

# handle the QUIETLY and REQUIRED arguments and set CULA_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "CULA" DEFAULT_MSG CULA_INCLUDE_DIR CULA_LIBRARIES )

MARK_AS_ADVANCED( CULA_INCLUDE_DIR CULA_LIBRARIES )
