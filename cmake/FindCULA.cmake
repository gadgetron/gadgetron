# - Find CULA
# Find the native CULA includes and library
#
#   CULA_FOUND       - True if TinyXML found.
#   CULA_INCLUDE_DIR - where to find tinyxml.h, etc.
#   CULA_LIBRARIES   - List of libraries when using TinyXML.
#

IF( CULA_INCLUDE_DIR )
    # Already in cache, be silent
    SET( TinyXML_FIND_QUIETLY TRUE )
ENDIF( CULA_INCLUDE_DIR )

FIND_PATH( CULA_INCLUDE_DIR "cula.h"
           PATH_SUFFIXES "cula/include" )

FIND_LIBRARY( CULA_LIBRARIES
              NAMES "cula"
              PATH_SUFFIXES "cula/lib" "cula/lib64" )

# handle the QUIETLY and REQUIRED arguments and set CULA_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "CULA" DEFAULT_MSG CULA_INCLUDE_DIR CULA_LIBRARIES )

MARK_AS_ADVANCED( CULA_INCLUDE_DIR CULA_LIBRARIES )
