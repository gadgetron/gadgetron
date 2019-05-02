if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUFFD__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_BINARY_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
  ${Boost_INCLUDE_DIR}
  ${ISMRMRD_INCLUDE_DIR}
  )

set(ffd_header_bspline_files 
                    BSplineFFD.h
                    BSplineFFD2D.h
                    BSplineFFD3D.h
                    BSplineFFD4D.h )

set(ffd_header_files ${ffd_header_bspline_files}
                    FFDBase.h
                    MLFFD.h)

set(ffd_src_files )

source_group(Base FILES FFDBase.h)
source_group(BSpline FILES ${ffd_header_bspline_files})
source_group(MLFFD FILES MLFFD.h)

#add_library(gadgetron_toolbox_cpuffd SHARED ${ffd_header_files} ${ffd_src_files} )

#set_target_properties(gadgetron_toolbox_cpuffd PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

#install(TARGETS gadgetron_toolbox_cpuffd DESTINATION lib COMPONENT main)

install(FILES
        ${ffd_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
