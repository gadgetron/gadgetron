if (CUDA_FOUND)
add_library(pingvin_toolbox_deblur SHARED
        MFIOperator.h
        MFIOperator.cpp)

set_target_properties(pingvin_toolbox_deblur PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})

target_link_libraries(pingvin_toolbox_deblur
        pingvin_toolbox_cpucore
        pingvin_toolbox_cpucore_math
        pingvin_toolbox_gpu
        )
target_include_directories(pingvin_toolbox_deblur
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

install(TARGETS pingvin_toolbox_deblur
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        MFIOperator.h
        DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)
endif()