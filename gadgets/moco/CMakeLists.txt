if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_MOCO__)
endif ()

if (CUDA_FOUND)
    list(APPEND CPU_GADGETS cpuRegistrationAveragingGadget.cpp)
    list(APPEND CPU_LIBS gadgetron_toolbox_cpucore_math gadgetron_toolbox_cpureg)


    list(APPEND GPU_GADGETS gpuRegistrationAveragingGadget.cpp gpuRegistrationScatteringGadget.cpp)
    list(APPEND GPU_LIBS gadgetron_toolbox_gpu ${CUDA_LIBRARIES})


    add_library(gadgetron_moco SHARED
            cpuRegistrationAveragingGadget.h
            gadgetron_moco_export.h
            gpuRegistrationAveragingGadget.h
            gpuRegistrationScatteringGadget.h
            RegistrationAveragingGadget.h
            RegistrationScatteringGadget.h
            ${CPU_GADGETS}
            ${GPU_GADGETS}
            )

    set_target_properties(gadgetron_moco PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

    target_link_libraries(gadgetron_moco
            gadgetron_core
            gadgetron_toolbox_cpucore gadgetron_mricore
            gadgetron_toolbox_gpu
            )

	install(TARGETS gadgetron_moco
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)
    install(TARGETS gadgetron_moco DESTINATION lib COMPONENT main)

    install(FILES cpuRegistrationAveragingGadget.h
            gadgetron_moco_export.h
            gpuRegistrationAveragingGadget.h
            gpuRegistrationScatteringGadget.h
            RegistrationAveragingGadget.h
            RegistrationScatteringGadget.h
            DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

    add_subdirectory(config)
elseif (CUDA_FOUND)
    message("Cuda not found, not compiling gpu-based registration gadgets")
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)