if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_PR__)
endif ()

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${ARMADILLO_INCLUDE_DIRS} 
    ${ACE_INCLUDE_DIR} 
    ${ISMRMRD_INCLUDE_DIR}
    ${FFTW3_INCLUDE_DIR}
    ${MKL_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/klt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/dwt/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/ffd
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
    ${CMAKE_SOURCE_DIR}/toolboxes/operators
    ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/image_io
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/pattern_recognition
  )

set(pr_header_fiels pr_export.h 
                    pr_kmeans.h )

set(pr_src_fiels pr_kmeans.cpp )

add_library(gadgetron_toolbox_pr SHARED 
            ${pr_header_fiels} 
            ${pr_src_fiels} )

set_target_properties(gadgetron_toolbox_pr PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_pr
                    ${ISMRMRD_LIBRARIES}
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_log
                    gadgetron_toolbox_cpucore_math
                    gadgetron_toolbox_cpufft 
                    gadgetron_toolbox_cpuklt 
                    gadgetron_toolbox_image_analyze_io 
                    gadgetron_toolbox_cpudwt 
                    gadgetron_toolbox_cpu_image )

install(TARGETS gadgetron_toolbox_pr DESTINATION lib COMPONENT main)

install(FILES ${pr_header_fiels} 
      DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
