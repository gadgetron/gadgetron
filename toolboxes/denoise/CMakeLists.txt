
include_directories(
        ${Boost_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${EIGEN_INCLUDE_DIR}
        ${MKL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_BINARY_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
)

add_library(gadgetron_toolbox_denoise SHARED 
        denoise_export.h
        non_local_means.cpp
        non_local_means.h
        non_local_bayes.cpp non_local_bayes.h)

set_target_properties(gadgetron_toolbox_denoise PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_denoise
        gadgetron_toolbox_log
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        ${ARMADILLO_LIBRARIES}
        ${Boost_LIBRARIES}
        )

install(TARGETS gadgetron_toolbox_denoise DESTINATION lib COMPONENT main)

install(FILES
        denoise_export.h
        non_local_means.h
        non_local_bayes.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
