set(pingvin_toolbox_python_header_files python_converters.h
        python_hoNDArray_converter.h
        python_mrd_converters.h
        python_numpy_wrappers.h
        python_toolbox.h
        python_tuple_converter.h
        python_optional_converter.h
        python_vector_converter.h)

add_library(pingvin_toolbox_python SHARED
        python_toolbox.cpp
        ${pingvin_toolbox_python_header_files})

target_link_libraries(pingvin_toolbox_python
        pingvin_core
        pingvin_toolbox_mri_core
        pingvin_toolbox_cpucore
        Boost::system
        python
        pingvin_toolbox_log)

target_include_directories(pingvin_toolbox_python
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)
set_target_properties(pingvin_toolbox_python PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

install(TARGETS pingvin_toolbox_python
	EXPORT pingvin-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        ${pingvin_toolbox_python_header_files}
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)