set( cpu_solver_header_files
        cpusolver_export.h
        hoGdSolver.h
        hoCgPreconditioner.h
        hoCgSolver.h
        hoLsqrSolver.h
        hoGpBbSolver.h
        hoSbCgSolver.h
        hoSolverUtils.h
        curveFittingSolver.h
        HybridLM.h
        simplexLagariaSolver.h )

add_library(gadgetron_toolbox_cpu_solver INTERFACE)
target_link_libraries(gadgetron_toolbox_cpu_solver INTERFACE
                    gadgetron_toolbox_operator
                    gadgetron_toolbox_solvers
                    gadgetron_toolbox_cpucore
                    gadgetron_toolbox_cpucore_math )
target_include_directories(gadgetron_toolbox_cpu_solver
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS gadgetron_toolbox_cpu_solver EXPORT gadgetron-export DESTINATION lib COMPONENT main)

install(FILES ${cpu_solver_header_files} 
    DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)