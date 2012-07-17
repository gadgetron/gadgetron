# - Find CodeSource Xsd
# This module can be used to find Xsd and it's include path
# Variables:
#	XSD_EXECUTABLE
#	XSD_INCLUDE_DIR
#	XSD_FOUND

SET(XSD_FOUND FALSE)

if(WIN32)
	SET(__XSD_NAME xsd.exe)
else(WIN32)
	SET(__XSD_NAME xsdcxx xsd)
endif(WIN32)

if(XSD_INCLUDE_DIR)
	#in cache already
	SET(XSD_FOUND TRUE)
else(XSD_INCLUDE_DIR)
	find_file(XSD_EXECUTABLE NAMES ${__XSD_NAME}
	     PATHS
		 ${XSD_DIR}/bin
		/usr/bin
		/usr/local/bin
	)

	if(XSD_EXECUTABLE)
		SET(XSD_FOUND TRUE)
	else(XSD_EXECUTABLE)
		SET(XSD_EXECUTABLE "xsd-NOTFOUND" CACHE FILE "xsd executable path")
	endif(XSD_EXECUTABLE)

	find_path(XSD_INCLUDE_DIR NAMES xsd 
	    PATHS
		${XSD_DIR}/include
		/usr/include
		/usr/local/include
	)

	if(XSD_INCLUDE_DIR)
		SET(XSD_FOUND TRUE)
	else(XSD_INCLUDE_DIR)
		SET(XSD_INCLUDE_DIR "xsd-include-NOTFOUND" CACHE PATH "xsd include path")
	endif(XSD_INCLUDE_DIR)
endif(XSD_INCLUDE_DIR)

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

