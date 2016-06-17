include_directories(
    ${CMAKE_SOURCE_DIR}/apps/gadgetron
    ${CMAKE_BINARY_DIR}/apps/gadgetron
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
    ${ISMRMRD_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${PYTHON_INCLUDE_PATH}
    ${NUMPY_INCLUDE_DIRS}
    ${ACE_INCLUDE_DIR})

add_library(gadgetron_toolbox_python SHARED
    python_toolbox.cpp
    python_toolbox.h
    python_export.h)

target_link_libraries(gadgetron_toolbox_python
    gadgetron_toolbox_cpucore
    ${ISMRMRD_LIBRARIES}
    ${PYTHON_LIBRARIES}
    ${Boost_LIBRARIES}
    ${ACE_LIBRARIES}
    gadgetron_toolbox_log)

set_target_properties(gadgetron_toolbox_python PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
if (WIN32)
    # only add this definition to the library, not the example binaries
    set_target_properties(gadgetron_toolbox_python PROPERTIES COMPILE_DEFINITIONS "__BUILD_GADGETRON_TOOLBOX_PYTHON__")
endif (WIN32)

install(TARGETS gadgetron_toolbox_python DESTINATION lib COMPONENT main)

install(FILES
    python_toolbox.h
    python_converters.h
    python_tuple_converter.h
    python_hoNDArray_converter.h
    python_ismrmrd_converter.h
    python_numpy_wrappers.h
    python_export.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

add_subdirectory(example)
