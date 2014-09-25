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
set(DEBIAN_PACKAGE_DEPENDS "libhdf5-serial-dev, git-core, libboost-all-dev, build-essential, libfftw3-dev, h5utils, hdf5-tools, libqt4-dev, libglew1.6-dev, xsltproc, fop, python-dev, python-numpy, freeglut3-dev, libxi-dev, liblapack-dev, libxml2-dev, libxslt-dev, libarmadillo-dev, libace-dev, python-h5py, python-matplotlib, python-libxml2, gcc-multilib, python-psutil, libgtest-dev, libboost-system-dev, libboost-thread-dev, libboost-program-options-dev, python-twisted")

# where the package metadata are
set(GADGETRON_CPACK_CFG_FILE "${PROJECT_BINARY_DIR}/cpack_options.cmake")

# where the package to be installed
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})