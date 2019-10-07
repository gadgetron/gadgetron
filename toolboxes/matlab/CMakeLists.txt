
add_library(gadgetron_toolbox_matlab SHARED  matlab_export.h MatlabUtils.h MatlabUtils.cpp )

set_target_properties(gadgetron_toolbox_matlab PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
target_include_directories(gadgetron_toolbox_matlab
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${MATLAB_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_core
        )


target_link_libraries(gadgetron_toolbox_matlab
        gadgetron_toolbox_mri_core
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        ${MATLAB_LIBRARIES}
        ${ISMRMRD_LIBRARIES})


install(TARGETS gadgetron_toolbox_matlab
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES MatlabUtils.h DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
