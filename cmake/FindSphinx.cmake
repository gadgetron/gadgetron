#Look for an executable called sphinx-build
find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build
             DOC "Path to sphinx-build executable")

             
#Look for the version of sphinx
execute_process(COMMAND ${SPHINX_EXECUTABLE} --version
                OUTPUT_VARIABLE _sphinx_version_output_string
                OUTPUT_STRIP_TRAILING_WHITESPACE)
#string(FIND ${_sphinx_version_output_string} " " _last_space_position REVERSE)
#math(EXPR _version_substring_begin "${_last_space_position} + 1")
#string(SUBSTRING ${_sphinx_version_output_string} ${_version_substring_begin} -1 SPHINX_VERSION)

include(FindPackageHandleStandardArgs)

#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args(Sphinx
                                  "Failed to find sphinx-build executable"
                                  SPHINX_EXECUTABLE)