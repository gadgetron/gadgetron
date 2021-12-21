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

# Enable/Disable automatic search for dependencies:
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# Enable/Disable component install for CPack generator DEB
set(CPACK_DEB_COMPONENT_INSTALL OFF)
set(CPACK_DEB_PACKAGE_COMPONENT OFF)

# Set dependencies explicitly
set(CPACK_DEBIAN_PACKAGE_DEPENDS "build-essential, ismrmrd, libfftw3-dev, python-dev, python-numpy, python-psutil, liblapack-dev, libxml2-dev, libxslt-dev, libarmadillo-dev, libace-dev, python-matplotlib, python-libxml2, python-h5py, libboost-all-dev, libhdf5-serial-dev, h5utils, hdf5-tools, libgtest-dev")

# Where the package metadata are
set(GADGETRON_CPACK_CFG_FILE "${PROJECT_BINARY_DIR}/cpack_options.cmake")

# Where the package to be installed
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
message("CPACK_PACKAGING_INSTALL_PREFIX: " ${CPACK_PACKAGING_INSTALL_PREFIX})

if(WIN32)
  # NSLS
  list(APPEND CPACK_GENERATOR "NSIS")    
endif()

list(APPEND CPACK_SOURCE_GENERATOR "TGZ")
list(APPEND CPACK_SOURCE_GENERATOR "ZIP")
list(APPEND CPACK_SOURCE_IGNORE_FILES ";.git;.gitignore;todo.txt;_clang-format;build/")

# Set dependencies explicitly
include(InstallRequiredSystemLibraries)
