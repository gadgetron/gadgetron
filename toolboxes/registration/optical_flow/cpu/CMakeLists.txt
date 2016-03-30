    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
    endif (WIN32)

    if(WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif(WIN32)

    include_directories(
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/transformation
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/solver
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/warper
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/dissimilarity
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/register
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu/application
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/algorithm
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
        ${CMAKE_SOURCE_DIR}/toolboxes/core/image_io
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
        ${CMAKE_SOURCE_DIR}/toolboxes/mri_image
        ${Boost_INCLUDE_DIR}
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

    set(transformation_files transformation/hoImageRegTransformation.h
                             transformation/hoImageRegParametricTransformation.h 
                             transformation/hoImageRegHomogenousTransformation.h 
                             transformation/hoImageRegRigid2DTransformation.h 
                             transformation/hoImageRegRigid3DTransformation.h 
                             transformation/hoImageRegNonParametricTransformation.h 
                             transformation/hoImageRegDeformationField.h )

    set(solver_files solver/hoImageRegSolver.h 
                     solver/hoImageRegParametricSolver.h
                     solver/hoImageRegParametricDownHillSolver.h
                     solver/hoImageRegParametricGradientDescentSolver.h
                     solver/hoImageRegNonParametricSolver.h
                     solver/hoImageRegDeformationFieldSolver.h 
                     solver/hoImageRegDeformationFieldBidirectionalSolver.h )

    set(warper_files warper/hoImageRegWarper.h)

    set(similarity_files dissimilarity/hoImageRegDissimilarity.h
                         dissimilarity/hoImageRegDissimilarityHistogramBased.h
                         dissimilarity/hoImageRegDissimilarityLocalCCR.h
                         dissimilarity/hoImageRegDissimilarityMutualInformation.h
                         dissimilarity/hoImageRegDissimilarityNormalizedMutualInformation.h
                         dissimilarity/hoImageRegDissimilaritySSD.h )

    set(register_files register/hoImageRegRegister.h
                       register/hoImageRegParametricRegister.h
                       register/hoImageRegNonParametricRegister.h
                       register/hoImageRegDeformationFieldRegister.h 
                       register/hoImageRegDeformationFieldBidirectionalRegister.h )

    set(application_files application/hoImageRegContainer2DRegistration.h )

    if ( BUILD_CPU_OPTIMAL_FLOW_REG )

        add_library(gadgetron_toolbox_cpureg SHARED 
                    cpureg_export.h
                    ${opticalflow_files} 
                    ${transformation_files} 
                    ${solver_files} 
                    ${warper_files}
                    ${similarity_files} 
                    ${register_files}
                    ${application_files} )

        set_target_properties(gadgetron_toolbox_cpureg PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})


        target_link_libraries(gadgetron_toolbox_cpureg 
                gadgetron_toolbox_cpucore
		gadgetron_toolbox_log
                gadgetron_toolbox_cpucore_math 
                ${ARMADILLO_LIBRARIES}
                optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} )

    else ( BUILD_CPU_OPTIMAL_FLOW_REG )
        message("cpureg library will not be built ... ")
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )

    if ( BUILD_CPU_OPTIMAL_FLOW_REG )
        source_group(opticalflow FILES ${opticalflow_files})
        install(FILES ${opticalflow_files} DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )

    source_group(transformation FILES ${transformation_files})
    source_group(solver FILES ${solver_files})
    source_group(warper FILES ${warper_files})
    source_group(similarity FILES ${similarity_files})
    source_group(register FILES ${register_files})
    source_group(application FILES ${application_files})

    if ( BUILD_CPU_OPTIMAL_FLOW_REG )
        install(TARGETS gadgetron_toolbox_cpureg DESTINATION lib COMPONENT main)
    endif ( BUILD_CPU_OPTIMAL_FLOW_REG )

    install(FILES
        ${transformation_files}
        ${solver_files}
        ${warper_files}
        ${similarity_files}
        ${register_files}
        ${application_files}
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
