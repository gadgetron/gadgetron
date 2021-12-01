
add_compile_options(-Werror)

if (CUDA_FOUND)
add_library(gadgetron_toolbox_deblur SHARED
        MFIOperator.h
        MFIOperator.cpp)

set_target_properties(gadgetron_toolbox_deblur PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_deblur
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_gpucore
        gadgetron_toolbox_gpuparallelmri
        gadgetron_toolbox_gpunfft
        )
target_include_directories(gadgetron_toolbox_deblur
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

install(TARGETS gadgetron_toolbox_deblur
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        MFIOperator.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
endif()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)