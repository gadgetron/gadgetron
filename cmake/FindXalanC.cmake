# - Try to find XalanC
# Once done this will define
#
#  XALANC_FOUND - System has XalanC
#  XALANC_INCLUDE_DIR - The XalanC include directory
#  XALANC_LIBRARY_DIR - The XalanC library dir
#  XALANC_LIBRARIES - The libraries needed to use XalanC
#  XALANC_DEFINITIONS - Compiler switches required for using XalanC

# Copyright (c) 2009, Helio Chissini de Castro, <helio@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.


IF (XALANC_INCLUDE_DIR AND XALANC_LIBRARIES)
   # in cache already
   SET(XalanC_FIND_QUIETLY TRUE)
ENDIF (XALANC_INCLUDE_DIR AND XALANC_LIBRARIES)


FIND_PATH(XALANC_INCLUDE_DIR DOMSupport/DOMServices.hpp
	PATHS
	/usr/local/include/xalanc
	/usr/include/xalanc
	PATH_SUFFIXES
	xalanc
	)

FIND_LIBRARY(XALANC_LIBRARIES NAMES xalan-c xalanMsg)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XalanC DEFAULT_MSG XALANC_LIBRARIES XALANC_INCLUDE_DIR)

MARK_AS_ADVANCED(XALANC_INCLUDE_DIR XALANC_LIBRARIES XALANC_LIBRARY_DIR)
