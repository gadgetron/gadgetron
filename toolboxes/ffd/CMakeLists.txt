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


add_library(pingvin_toolbox_cpuffd INTERFACE)

target_include_directories(pingvin_toolbox_cpuffd
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})
#
#add_library(pingvin_toolbox_cpuffd SHARED ${ffd_header_files} ${ffd_src_files} )
#
#set_target_properties(pingvin_toolbox_cpuffd PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})
#
#target_link_libraries(pingvin_toolbox_cpuffd
#        pingvin_toolbox_cpucore
#        )

#target_include_directories(
#        pingvin_toolbox_cpuffd
#        PUBLIC
#        ${CMAKE_CURRENT_SOURCE_DIR}
#)
#install(TARGETS pingvin_toolbox_cpuffd DESTINATION lib COMPONENT main)

install(FILES
        ${ffd_header_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)