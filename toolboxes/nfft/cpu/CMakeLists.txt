if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUNFFT__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

add_library(gadgetron_toolbox_cpunfft SHARED
    hoNFFT.h
    hoNFFT.cpp
    hoNFFT_sparseMatrix.h
    hoNFFT_sparseMatrix.cpp
    hoNFFTOperator.cpp
  )

set_target_properties(gadgetron_toolbox_cpunfft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_include_directories(gadgetron_toolbox_cpunfft
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/.. )
target_link_libraries(gadgetron_toolbox_cpunfft
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_log
    gadgetron_toolbox_hostutils
    gadgetron_toolbox_operator
  )

install(TARGETS gadgetron_toolbox_cpunfft
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)

install(FILES 
    hoNFFT.h hoNFFT_sparseMatrix.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)