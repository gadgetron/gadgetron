if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CPUNFFT__)
  add_definitions(-D_USE_MATH_DEFINES)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()

include_directories( 
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${Boost_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
  )

add_library(gadgetron_toolbox_cpunfft SHARED 
    hoNFFT.h
    cpunfft_export.h
    hoNFFT.cpp
  )

set_target_properties(gadgetron_toolbox_cpunfft PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_cpunfft 
    gadgetron_toolbox_cpufft
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_log
    gadgetron_toolbox_hostutils
  )

install(TARGETS gadgetron_toolbox_cpunfft DESTINATION lib COMPONENT main)

install(FILES 
    hoNFFT.h 
    cpunfft_export.h 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
