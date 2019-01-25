
include_directories(
        ${Boost_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_BINARY_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
)

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
        ${ARMADILLO_LIBRARIES}
        ${Boost_LIBRARIES}
        )

install(TARGETS gadgetron_toolbox_deblur DESTINATION lib COMPONENT main)

install(FILES
        MFIOperator.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
