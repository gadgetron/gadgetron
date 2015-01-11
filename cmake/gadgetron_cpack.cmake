################################################################################
# Find available package generators
################################################################################

if(UNIX)
  # DEB
  find_program(DPKG_PROGRAM dpkg)
  if(EXISTS ${DPKG_PROGRAM})
    list(APPEND CPACK_GENERATOR "DEB")
  endif(EXISTS ${DPKG_PROGRAM})
endif(UNIX)

if(WIN32)
    # NSLS
    list(APPEND CPACK_GENERATOR "NSIS")    
endif(WIN32)

list(APPEND CPACK_SOURCE_GENERATOR "TGZ")
list(APPEND CPACK_SOURCE_GENERATOR "ZIP")
list(APPEND CPACK_SOURCE_IGNORE_FILES ";.git;.gitignore;todo.txt;_clang-format;build/")

# set dependencies explictly
include(InstallRequiredSystemLibraries)

# autogenerate dependency information
#set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(DEBIAN_PACKAGE_DEPENDS "libfftw3-dev, python, python-numpy, liblapack-dev, libxml2-dev, libxslt-dev, libarmadillo-dev, libace-dev, python-matplotlib, python-libxml2, libboost-system-dev, libboost-thread-dev, libboost-program-options-dev, libboost-chrono-dev, libboost-filesystem-dev, ismrmrd")

set(CPACK_COMPONENT_TESTCLIENT_DEPENDS "ismrmrd, libboost-program-options-dev, libboost-thread-dev, libboost-system-dev")

#set(CPACK_COMPONENT_MAIN_DEPENDS "libfftw3-dev, python, python-numpy, liblapack-dev, libxml2-dev, libxslt-dev, libarmadillo-dev, libace-dev, python-matplotlib, python-libxml2, libboost-system-dev, libboost-thread-dev, libboost-program-options-dev, libboost-chrono-dev, libboost-filesystem-dev, ismrmrd")

#set(CPACK_COMPONENT_TESTCLIENT_DEPENDS "main")

# where the package metadata are
set(GADGETRON_CPACK_CFG_FILE "${PROJECT_BINARY_DIR}/cpack_options.cmake")

# where the package to be installed
# set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX})
if (NOT WIN32)
    set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif (NOT WIN32)
