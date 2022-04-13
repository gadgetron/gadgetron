if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUGADGET__)
endif ()

add_library(gadgetron_gpugadget SHARED
        cuFFTGadget.h
        cuFFTGadget.cpp)

set_target_properties(gadgetron_gpugadget PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_gpugadget
        gadgetron_core
        gadgetron_toolbox_gpu)

install(FILES
        cuFFTGadget.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install(TARGETS gadgetron_gpugadget
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)