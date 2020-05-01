if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUFFT__)
endif ()


add_library(gadgetron_toolbox_cpufft SHARED
        cpufft_export.h
        hoNDFFT.h
        hoNDFFT.cpp
        )

set_target_properties(gadgetron_toolbox_cpufft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(
        gadgetron_toolbox_cpufft
        PUBLIC
        $<INSTALL_INTERFACE:include/gadgetron>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(gadgetron_toolbox_cpufft
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore_math
        FFTW
        armadillo
        boost
        )



install(TARGETS gadgetron_toolbox_cpufft
		EXPORT gadgetron-export
		LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES
        cpufft_export.h
        hoNDFFT.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)