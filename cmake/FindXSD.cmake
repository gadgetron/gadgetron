# - Find CodeSynthesis XSD
# This module can be used to find XSD and it's include path
# Variables:
#	XSD_FOUND - System has XSD
#	XSD_EXECUTABLE - XSD binary executable
#	XSD_INCLUDE_DIR - XSD include directory
#
# Functions:
#       WRAP_XSD - Generates C++ bindings in the given output directory for a given schema file

if(NOT DEFINED XSD_DIR AND DEFINED ENV{XSD_DIR})
    set(XSD_DIR $ENV{XSD_DIR})
endif(NOT DEFINED XSD_DIR AND DEFINED ENV{XSD_DIR})

find_program(XSD_EXECUTABLE NAMES xsd xsdcxx xsd.exe
    PATHS ${XSD_DIR} /usr /usr/local
    PATH_SUFFIXES bin
)

find_path(XSD_INCLUDE_DIR NAMES xsd/cxx/pre.hxx
    PATHS ${XSD_DIR} /usr /usr/local
    PATH_SUFFIXES include
)

FUNCTION(XSD_EXTRACT_OPTIONS _xsd_files _xsd_options)
	foreach(current_arg ${ARGN})
		IF(${current_arg} STREQUAL "OPTIONS")
			SET(_XSD_DOING_OPTIONS TRUE)
		else(${current_arg} STREQUAL "OPTIONS")
			if(_XSD_DOING_OPTIONS)
				SET(_xsd_options_p ${_xsd_options_p} ${current_arg})
			else(_XSD_DOING_OPTIONS)
				SET(_xsd_files_p ${_xsd_files_p} ${current_arg})
			endif(_XSD_DOING_OPTIONS)
		endif(${current_arg} STREQUAL "OPTIONS")
	endforeach(current_arg)
	SET(${_xsd_files} ${_xsd_files_p} PARENT_SCOPE)
	SET(${_xsd_options} ${_xsd_options_p} PARENT_SCOPE)
ENDFUNCTION(XSD_EXTRACT_OPTIONS)


FUNCTION(WRAP_XSD XSD_SRCS XSD_INCLUDES OUT_PATH)
	SET(OUTPUT_DIR  ${CMAKE_CURRENT_BINARY_DIR}/src/xsd)
	FILE(MAKE_DIRECTORY ${OUTPUT_DIR})
	SET(${XSD_INCLUDES} ${OUTPUT_DIR} PARENT_SCOPE)
	XSD_EXTRACT_OPTIONS(xsd_files xsd_options ${ARGN})
	FOREACH(it ${xsd_files})
		STRING(REGEX REPLACE ".*/" "" BARE_XSD "${it}" )
		STRING(REGEX REPLACE ".xsd" ".cxx" SOURCE "${BARE_XSD}" )
		STRING(REGEX REPLACE ".xsd" ".hxx" HEADER "${BARE_XSD}" )
		CONFIGURE_FILE(${it} ${OUT_PATH}/${BARE_XSD} COPY_ONLY)
		SET(SOURCE ${OUTPUT_DIR}/${SOURCE})
		SET(HEADER ${OUTPUT_DIR}/${HEADER})
		ADD_CUSTOM_COMMAND(OUTPUT ${SOURCE} ${HEADER}
				COMMAND ${XSD_EXECUTABLE} ${xsd_options} "--output-dir" ${OUTPUT_DIR} ${OUT_PATH}/${BARE_XSD}
				DEPENDS ${it}
				VERBATIM
		)
		set_source_files_properties(${HEADER} PROPERTIES GENERATED TRUE)
		set_source_files_properties(${SOURCE} PROPERTIES GENERATED TRUE)
		SET(_XSD_SRCS ${_XSD_SRCS} ${SOURCE} ${HEADER})
	ENDFOREACH(it)
	SET(${XSD_SRCS} ${_XSD_SRCS} PARENT_SCOPE)
ENDFUNCTION(WRAP_XSD)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XSD DEFAULT_MSG XSD_INCLUDE_DIR XSD_EXECUTABLE)
mark_as_advanced(XSD_INCLUDE_DIR XSD_EXECUTABLE)
