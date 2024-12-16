find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
    find_package(Sphinx REQUIRED)

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)
    add_custom_target(Sphinx ALL
            COMMAND ${SPHINX_EXECUTABLE} -b html
            # Tell Breathe where to find the Doxygen output
            -Dbreathe_projects.Pingvin=${CMAKE_CURRENT_BINARY_DIR}/doxygen/doc/xml/
            ${SPHINX_SOURCE} ${SPHINX_BUILD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating documentation with Sphinx")
    add_dependencies(Sphinx Doxygen)
else()
	message("Doxygen not found. Will not be able to build documentation")
endif()
