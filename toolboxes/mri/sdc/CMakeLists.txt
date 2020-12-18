add_subdirectory(cpu)

if (CUDA_FOUND)
  add_subdirectory(gpu)
endif ()

add_library(gadgetron_toolbox_sdc INTERFACE)
target_include_directories(gadgetron_toolbox_sdc
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR})

install(FILES
    SDC.h
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)