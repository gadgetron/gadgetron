find_package(Doxygen)
if(DOXYGEN_FOUND)
	add_subdirectory(doxygen)
else(DOXYGEN_FOUND)
	MESSAGE("Doxygen not found. Will not be able to build documentation")
endif(DOXYGEN_FOUND)

find_file(XSLTPROCEXE xsltproc
		HINTS /usr/bin
		/usr/local/bin)

find_path(DOCBOOK_XSL_DIR html/docbook.xsl
  PATHS /usr/local/share/xml/xls/docbook-xsl-ns
  		/usr/share/xml/xls/docbook-xsl-ns
  		/usr/share/xml/docbook 
 		/usr/share/xml/docbook/stylesheet/docbook-xsl-ns/
       $ENV{DOCBOOKDIR}
  NO_DEFAULT_PATH)

if (NOT DOCBOOK_XSL_DIR)
  message("Could not find HTML docbook.xsl, try to set DOCBOOKDIR")
endif (NOT DOCBOOK_XSL_DIR)

if (XSLTPROCEXE AND DOCBOOK_XSL_DIR)
	MESSAGE("Docbook tools found, building manual XSLTPROCEXE: ${XSLTPROCEXE}")
	add_subdirectory(manual)
else (XSLTPROCEXE AND DOCBOOK_XSL_DIR)
    MESSAGE("Docbook tools not found")
endif (XSLTPROCEXE AND DOCBOOK_XSL_DIR)