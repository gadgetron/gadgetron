if (WIN32)
    add_definitions(-D__BUILD_GADGETRON_CPUSOLVERS__)
endif ()

include_directories(
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
    ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers 
    ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu )

set( cpu_solver_header_files
        ../cgCallback.h
        ../cgPreconditioner.h
        ../cgSolver.h
        ../eigenTester.h
        ../gpBbSolver.h
        ../gpSolver.h
        ../lbfgsSolver.h
        ../linearOperatorSolver.h
        ../lsqrSolver.h
        ../lwSolver.h
        ../nlcgSolver.h
        ../osLALMSolver.h
        ../osMOMSolver.h
        ../osSPSSolver.h
        ../sbcSolver.h
        ../sbSolver.h
        ../solver.h

        cpusolver_export.h
        hoGdSolver.h 
        hoCgPreconditioner.h 
        hoCgSolver.h 
        hoLsqrSolver.h 
        hoGpBbSolver.h 
        hoSbCgSolver.h 
        hoSolverUtils.h 
        curveFittingSolver.h 
        simplexLagariaSolver.h )

set( cpu_solver_source_files )

add_library(gadgetron_toolbox_cpu_solver SHARED ${cpu_solver_header_files} ${cpu_solver_source_files} )

set_target_properties(gadgetron_toolbox_cpu_solver PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})
set_target_properties(gadgetron_toolbox_cpu_solver PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(gadgetron_toolbox_cpu_solver 
                    gadgetron_toolbox_cpucore 
                    gadgetron_toolbox_cpucore_math )

# install(TARGETS gadgetron_toolbox_cpu_solver DESTINATION lib COMPONENT main)

install(FILES ${cpu_solver_header_files} 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
