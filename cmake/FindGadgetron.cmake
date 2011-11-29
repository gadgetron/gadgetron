#
# Find the Gadgetron Installation
#

# This module defines
# GADGETRON_INCLUDE_DIR, where to finds Gadget.h
# GADGETRON_HOME, Gadgetron Root Dir
# GADGETRON_LIB_DIR, This is where all the installed gadgetron libraries live
# GADGETRON_FOUND, if false, you cannot build anything that requires ACE

# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set(_check_list)

# Search for the header file.
find_path(GADGETRON_HOME include/Gadget.h 
  HINTS $ENV{GADGETRON_HOME} /usr/local/gadgetron /usr/gadgetron)
mark_as_advanced(GADGETRON_HOME)
list(APPEND _check_list GADGETRON_HOME)

SET(GADGETRON_INCLUDE_DIR ${GADGETRON_HOME}/include)
mark_as_advanced(GADGETRON_INCLUDE_DIR)
list(APPEND _check_list GADGETRON_INCLUDE_DIR)

SET(GADGETRON_LIB_DIR ${GADGETRON_HOME}/lib)
mark_as_advanced(GADGETRON_LIB_DIR)
list(APPEND _check_list GADGETRON_LIB_DIR)

# Handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gadgetron DEFAULT_MSG ${_check_list})