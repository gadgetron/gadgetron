find_library(PUGIXML_LIBRARIES
        NAMES pugixml
        )
find_path(PUGIXML_INCLUDE_DIRS
        NAMES pugixml.hpp
        PATH_SUFFIXES pugixml
        )

# Read and parse pugixml version header file for version number
file(STRINGS "${PUGIXML_INCLUDE_DIRS}/pugixml.hpp" _pugixml_HEADER_CONTENTS REGEX "PUGIXML_VERSION")
string(REGEX REPLACE ".*PUGIXML_VERSION (.*)" "\\1" PUGIXML_VERSION "${_pugixml_HEADER_CONTENTS}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PugiXML DEFAULT_MSG
        PUGIXML_LIBRARIES
        PUGIXML_INCLUDE_DIRS
        )