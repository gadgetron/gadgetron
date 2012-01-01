

add_custom_target(htmlmanual ${XSLTPROCEXE}
${DOCBOOK_XSL_DIR}/html/docbook.xsl ${CMAKE_CURRENT_SOURCE_DIR}/gadgetron_manual.xml > ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.html WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating HTML docbook manual" VERBATIM)

find_file(FOPEXE fop
		HINTS /usr/bin
		/usr/local/bin)

if (FOPEXE)
	MESSAGE("FOP found, Building PDF Manual, FOPEXE: ${FOPEXE}")
	add_custom_target(fomanual ${XSLTPROCEXE}
	${DOCBOOK_XSL_DIR}/fo/docbook.xsl ${CMAKE_CURRENT_SOURCE_DIR}/gadgetron_manual.xml > ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.fo WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generating FO (for PDF conversion) docbook manual" VERBATIM)
	
	add_custom_target(pdfmanual ${FOPEXE} ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.fo ${CMAKE_CURRENT_BINARY_DIR}/gadgetron_manual.pdf DEPENDS fomanual)
	
else (FOPEXE)
	MESSAGE("FOP executable not found, PDF manual cannot be build")
endif(FOPEXE)