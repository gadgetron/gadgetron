

set(gadgetron_toolbox_python_header_files python_converters.h
        python_export.h
        python_hoNDArray_converter.h
        python_IsmrmrdReconData_converter.h
        python_IsmrmrdImageArray_converter.h
        python_ismrmrd_converter.h
        python_numpy_wrappers.h
        python_toolbox.h
        python_tuple_converter.h
        python_vector_converter.h)

add_library(gadgetron_toolbox_python SHARED
        python_toolbox.cpp
        ${gadgetron_toolbox_python_header_files})

message("PYTHON_INCLUDE_PATH is ${PYTHON_INCLUDE_PATH}")
message("PYTHON_LIBRARIES is ${PYTHON_LIBRARIES}")
message("NUMPY_INCLUDE_DIRS is ${NUMPY_INCLUDE_DIRS}")
message("Boost_LIBRARIES is ${Boost_LIBRARIES}")

target_link_libraries(gadgetron_toolbox_python
        gadgetron_core
        gadgetron_toolbox_mri_core
        gadgetron_toolbox_cpucore
        Boost::system
        python
        gadgetron_toolbox_log)

target_include_directories(gadgetron_toolbox_python
        PUBLIC
        $<INSTALL_INTERFACE:include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
		)
set_target_properties(gadgetron_toolbox_python PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
if (WIN32)
    # only add this definition to the library, not the example binaries
    set_target_properties(gadgetron_toolbox_python PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_TOOLBOX_PYTHON__")
endif ()

install(TARGETS gadgetron_toolbox_python
	EXPORT gadgetron-export
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        ${gadgetron_toolbox_python_header_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(example)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)