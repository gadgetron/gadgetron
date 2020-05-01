if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUDWT__)
endif ()



add_library(gadgetron_toolbox_cpudwt SHARED 
    cpudwt_export.h 
    hoNDWavelet.h
    hoNDWavelet.cpp
    hoNDHarrWavelet.h
    hoNDHarrWavelet.cpp
    hoNDRedundantWavelet.h
    hoNDRedundantWavelet.cpp
    )
target_include_directories(gadgetron_toolbox_cpudwt
        PUBLIC

        $<INSTALL_INTERFACE:include/gadgetron>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

set_target_properties(gadgetron_toolbox_cpudwt PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cpudwt
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math
                    )

install(TARGETS gadgetron_toolbox_cpudwt
        EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
    cpudwt_export.h 
    hoNDWavelet.h
    hoNDHarrWavelet.h
    hoNDRedundantWavelet.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)