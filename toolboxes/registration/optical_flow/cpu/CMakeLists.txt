    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
    endif (WIN32)

    if(WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif(WIN32)

    include_directories(
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${CMAKE_SOURCE_DIR}/gadgets/core
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/util
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/workflow
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/gtplus/solver
        ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
        ${CMAKE_SOURCE_DIR}/apps/gadgetron
        ${CMAKE_SOURCE_DIR}/apps/matlab
        ${CMAKE_SOURCE_DIR}/gadgets/mri_core 
        ${CMAKE_SOURCE_DIR}/gadgets/gtPlus
        ${Boost_INCLUDE_DIR}
        ${MKL_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIR} 
        ${ACE_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR} )

    set(opticalflow_files 
                hoCKOpticalFlowSolver.cpp
                hoCKOpticalFlowSolver.h
                hoHSOpticalFlowSolver.cpp
                hoHSOpticalFlowSolver.h
                hoLinearResampleOperator.cpp
                hoLinearResampleOperator.h
                hoOpticalFlowSolver.cpp
                hoOpticalFlowSolver.h )

    if ( BUILD_CPU_OPTIMAL_FLOW_REG )

        add_library(cpureg ${LIBTYPE} 
                    cpureg_export.h
                    ${opticalflow_files} )

        target_link_libraries(cpureg 
                cpucore
                cpucore_math 
                gtplus
                ${ARMADILLO_LIBRARIES}
                optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

    else ( BUILD_CPU_OPTIMAL_FLOW_REG )
        message("cpureg library will not be built ... ")
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )
   
    if ( BUILD_CPU_OPTIMAL_FLOW_REG )
        source_group(opticalflow FILES ${opticalflow_files})
        install(FILES ${opticalflow_files} DESTINATION include)
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )

    if ( BUILD_CPU_OPTIMAL_FLOW_REG )
        install(TARGETS cpureg DESTINATION lib)
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )
