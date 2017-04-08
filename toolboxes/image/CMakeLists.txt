if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_IMAGE__)
endif ()

include_directories(
    ${Boost_INCLUDE_DIR} 
    ${CMAKE_SOURCE_DIR}/toolboxes/core
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
    ${CMAKE_SOURCE_DIR}/toolboxes/image
    ${CMAKE_SOURCE_DIR}/toolboxes/image/cpu
    )

add_subdirectory(cpu)