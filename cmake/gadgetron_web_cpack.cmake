################################################################################
# Find available package generators
################################################################################

if(UNIX)
  # DEB
  find_program(DPKG_PROGRAM dpkg)
  if(EXISTS ${DPKG_PROGRAM})
    list(APPEND CPACK_GENERATOR "DEB")
  endif()
endif()

if(WIN32)
    # NSLS
    list(APPEND CPACK_GENERATOR "NSIS")    
endif()

list(APPEND CPACK_SOURCE_GENERATOR "TGZ")
list(APPEND CPACK_SOURCE_GENERATOR "ZIP")
list(APPEND CPACK_SOURCE_IGNORE_FILES ";.git;.gitignore;todo.txt;_clang-format;build/")

# set dependencies explicitly
set(DEBIAN_PACKAGE_DEPENDS "gadgetron, python-psutil, python-twisted")

# where the package metadata are
set(GADGETRON_WEB_CPACK_CFG_FILE "${PROJECT_BINARY_DIR}/cpack_options_web.cmake")

# where the package to be installed
# set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX})
if (NOT WIN32)
    set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif ()
