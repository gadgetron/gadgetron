################################################################################
# Metadata for package generators
################################################################################

# Common options
set(CPACK_PACKAGE_VERSION "4.2.1")
set(CPACK_PACKAGE_VERSION_MAJOR "4")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_NAME "GADGETRON")
set(CPACK_PACKAGE_VENDOR "https://gadgetron.github.io/")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Gadgetron framework")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "gadgetron")
set(CPACK_RESOURCE_FILE_LICENSE "/home/javeda2/mrprogs/gadgetron_fork/LICENSE")
set(CPACK_PACKAGE_MAINTAINER "Michael S. Hansen <michael.hansen@nih.gov>")
set(CPACK_PACKAGE_CONTACT "Michael S. Hansen <michael.hansen@nih.gov>")

# DEB specific
set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Implementation of the Gadgetron.")
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "/home/javeda2/mrprogs/gadgetron_fork/cmake/debian/postinst;/home/javeda2/mrprogs/gadgetron_fork/cmake/debian/prerm;" )

# NSIS specific
set(CPACK_NSIS_HELP_LINK "https://github.com/gadgetron/gadgetron")
set(CPACK_NSIS_URL_INFO_ABOUT "https://github.com/gadgetron/gadgetron")
set(CPACK_NSIS_MODIFY_PATH ON)
set(CPACK_NSIS_DISPLAY_NAME "gadgetron")

set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "ExecWait '$INSTDIR/cmake/InstallWinGadgetron.bat'")

# Output filename of the generated tarball / package
set(CPACK_PACKAGE_FILE_NAME "gadgetron-4.2.1")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "gadgetron-4.2.1")
