
add_library(gadgetron_toolbox_demons demons_registration.cpp)
target_link_libraries(gadgetron_toolbox_demons gadgetron_toolbox_cpucore)


set_target_properties(gadgetron_toolbox_mri_core PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_mri_core PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(
        gadgetron_toolbox_demons
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)




install(FILES
        demons_registration.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

