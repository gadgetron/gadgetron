if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_GPUREG__)
endif ()


cuda_add_library(gadgetron_toolbox_gpureg SHARED
        cuOpticalFlowSolver.cu
        cuHSOpticalFlowSolver.cu
        cuCKOpticalFlowSolver.cu
        cuResampleOperator.cu
        cuLinearResampleOperator.cu
        )

set_target_properties(gadgetron_toolbox_gpureg PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpureg
        gadgetron_toolbox_registration_optical_flow
        gadgetron_toolbox_operator
        gadgetron_toolbox_gpucore
        gadgetron_toolbox_solvers
        gadgetron_toolbox_log
        ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

target_include_directories(gadgetron_toolbox_gpureg
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

install(TARGETS gadgetron_toolbox_gpureg
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	RUNTIME DESTINATION bin
	COMPONENT main
)


install(FILES
        cuOpticalFlowSolver.h
        cuHSOpticalFlowSolver.h
        cuCKOpticalFlowSolver.h
        gpureg_export.h
        cuResampleOperator.h
        cuLinearResampleOperator.h
        cuCGHSOFSolver.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)