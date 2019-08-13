if (DOXYGEN_FOUND)
    set(DOXYGEN_FILE ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${DOXYGEN_FILE} @ONLY)
    add_custom_target(apidoc ${DOXYGEN_EXECUTABLE}
            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM)


else ()
    message("Doxygen not found. Will not be able to build documentation")
endif ()