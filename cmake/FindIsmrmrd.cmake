# - Find ISMRMRRD
#   ISMRMRD_FOUND            - true if an ISMRMRD installation is found.
#   ISMRMRD_INCLUDE_DIR      - where to find ismrmrd.h, etc.
#   ISMRMRD_LIBRARIES        - libismrmrd.so and libismrmrd_xml.so
#   ISMRMRD_SCHEMA_DIR       - where to find ismrmrd.xsd       

FIND_PATH( ISMRMRD_INCLUDE_DIR ismrmrd/ismrmrd.h 
HINTS $ENV{ISMRMRD_HOME} PATHS /usr/local /usr PATH_SUFFIXES include)

FIND_PATH( ISMRMRD_SCHEMA_DIR ismrmrd.xsd 
HINTS $ENV{ISMRMRD_HOME} PATHS /usr/local /usr PATH_SUFFIXES share/ismrmrd/schema)

FIND_LIBRARY( ISMRMRD_LIBRARY NAMES ismrmrd
HINTS $ENV{ISMRMRD_HOME} /usr/local /usr PATH_SUFFIXES lib)

SET(ISMRMRD_LIBRARIES ${ISMRMRD_LIBRARY})

INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "Ismrmrd" DEFAULT_MSG ISMRMRD_INCLUDE_DIR ISMRMRD_LIBRARIES ISMRMRD_SCHEMA_DIR)

MARK_AS_ADVANCED( ISMRMRD_INCLUDE_DIR ISMRMRD_LIBRARIES ISMRMRD_SCHEMA_DIR)

#if(ISMRMRD_FOUND)
#  message("ISMRMRD found ${ISMRMRD_LIBRARIES}")
#endif(ISMRMRD_FOUND)
