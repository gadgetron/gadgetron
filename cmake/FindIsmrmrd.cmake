# - Find ISMRMRRD
#   ISMRMRD_FOUND        - True if ISMRMRD found.
#   ISMRMRD_INCLUDE_DIR  - where to find ismrmrd.h, etc.
#   ISMRMRD_LIBRARIES    - libismrmrd.so, etc.
#   ISMRMRD_SCHEMA_DIR   - Where to find ismrmrd.xsd       

FIND_PATH( ISMRMRD_INCLUDE_DIR ismrmrd.h PATHS /usr/local/ /usr/include /usr/local/include PATH_SUFFIXES ismrmrd ismrmrd/include)
FIND_PATH( ISMRMRD_SCHEMA_DIR  ismrmrd.xsd PATHS /usr/local/ /usr/include /usr/local/include PATH_SUFFIXES ismrmrd ismrmrd/schema)

FIND_LIBRARY( ISMRMRD_LIBRARIES
              NAMES "ismrmrd"
              PATHS  /usr/local/lib ${ISMRMRD_INCLUDE_DIR}/../lib /usr/lib )

INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "Ismrmrd" DEFAULT_MSG ISMRMRD_INCLUDE_DIR ISMRMRD_LIBRARIES ISMRMRD_SCHEMA_DIR)

MARK_AS_ADVANCED( ISMRMRD_INCLUDE_DIR ISMRMRD_LIBRARIES ISMRMRD_SCHEMA_DIR)
