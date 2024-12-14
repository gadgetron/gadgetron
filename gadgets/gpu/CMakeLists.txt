add_library(pingvin_gpugadget SHARED
        cuFFTGadget.h
        cuFFTGadget.cpp)

set_target_properties(pingvin_gpugadget PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_gpugadget
        pingvin_core
        pingvin_toolbox_gpu)

install(FILES
        cuFFTGadget.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS pingvin_gpugadget
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)