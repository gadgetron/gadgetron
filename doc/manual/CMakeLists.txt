

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/figs ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)

#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/figs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_custom_target(htmlmanual ${XSLTPROCEXE} --xinclude -o ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.html
${DOCBOOK_XSL_DIR}/html/docbook.xsl ${CMAKE_CURRENT_SOURCE_DIR}/gadgetron_manual.xml  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating HTML docbook manual" VERBATIM)

add_custom_command(TARGET htmlmanual POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/figs ${CMAKE_CURRENT_BINARY_DIR}/figs
	COMMENT  "Copying figure files" VERBATIM)

find_file(FOPEXE fop
		HINTS /usr/bin
		/usr/local/fop
		/usr/local/bin)

if (FOPEXE)
	MESSAGE("FOP found, Building PDF Manual, FOPEXE: ${FOPEXE}")
	
	add_custom_target(pdfmanual ${FOPEXE} ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.fo ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.pdf  DEPENDS fomanual htmlmanual
	COMMENT "Generating PDF docbook manual" VERBATIM)
	
	add_custom_command(OUTPUT fomanual COMMAND ${XSLTPROCEXE} --xinclude -o ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.fo
	${DOCBOOK_XSL_DIR}/fo/docbook.xsl ${CMAKE_CURRENT_SOURCE_DIR}/gadgetron_manual.xml  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generating FO file for PDF docbook manual" VERBATIM)
	
	
else (FOPEXE)
	MESSAGE("FOP executable not found, PDF manual cannot be build")
endif(FOPEXE)