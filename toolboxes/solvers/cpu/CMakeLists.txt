if (WIN32)
    ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUSOLVERS__)
endif (WIN32)

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu )

set( cpu_solver_header_files
        cpusolver_export.h
        hoCgPreconditioner.h 
        hoCgSolver.h 
        hoLSQRSolver.h 
        hoGpBbSolver.h 
        hoSbCgSolver.h 
        hoSolverUtils.h )

set( cpu_solver_source_files
        hoLSQRSolver.cpp )

add_library(gadgetron_toolbox_cpu_solver SHARED ${cpu_solver_header_files} ${cpu_solver_source_files} )

set_target_properties(gadgetron_toolbox_cpu_solver PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_cpu_solver PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_cpu_solver 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math )

install(TARGETS gadgetron_toolbox_cpu_solver DESTINATION lib COMPONENT main)

install(FILES ${cpu_solver_header_files} 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
