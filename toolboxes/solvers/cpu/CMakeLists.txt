set( cpu_solver_header_files
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

add_library(pingvin_toolbox_cpu_solver INTERFACE)
target_link_libraries(pingvin_toolbox_cpu_solver INTERFACE
                    pingvin_toolbox_operator
                    pingvin_toolbox_solvers
                    pingvin_toolbox_cpucore
                    pingvin_toolbox_cpucore_math )
target_include_directories(pingvin_toolbox_cpu_solver
        INTERFACE
        $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS pingvin_toolbox_cpu_solver EXPORT pingvin-export DESTINATION lib COMPONENT main)

install(FILES ${cpu_solver_header_files}
    DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)