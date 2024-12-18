    set(opticalflow_files
            hoCKOpticalFlowSolver.cpp
            hoCKOpticalFlowSolver.h
            hoHSOpticalFlowSolver.cpp
            hoHSOpticalFlowSolver.h
            hoOpticalFlowSolver.cpp
            hoOpticalFlowSolver.h
            hoLinearResampleOperator.cpp
            hoLinearResampleOperator.h
            )

    set(transformation_files transformation/hoImageRegTransformation.h
            transformation/hoImageRegParametricTransformation.h
            transformation/hoImageRegHomogenousTransformation.h
            transformation/hoImageRegRigid2DTransformation.h
            transformation/hoImageRegRigid3DTransformation.h
            transformation/hoImageRegNonParametricTransformation.h
            transformation/hoImageRegDeformationField.h)

    set(solver_files solver/hoImageRegSolver.h
            solver/hoImageRegParametricSolver.h
            solver/hoImageRegParametricDownHillSolver.h
            solver/hoImageRegParametricGradientDescentSolver.h
            solver/hoImageRegNonParametricSolver.h
            solver/hoImageRegDeformationFieldSolver.h
            solver/hoImageRegDeformationFieldBidirectionalSolver.h)

    set(warper_files warper/hoImageRegWarper.h)

    set(similarity_files dissimilarity/hoImageRegDissimilarity.h
            dissimilarity/hoImageRegDissimilarityHistogramBased.h
            dissimilarity/hoImageRegDissimilarityLocalCCR.h
            dissimilarity/hoImageRegDissimilarityMutualInformation.h
            dissimilarity/hoImageRegDissimilarityNormalizedMutualInformation.h
            dissimilarity/hoImageRegDissimilaritySSD.h)

    set(register_files register/hoImageRegRegister.h
            register/hoImageRegParametricRegister.h
            register/hoImageRegNonParametricRegister.h
            register/hoImageRegDeformationFieldRegister.h
            register/hoImageRegDeformationFieldBidirectionalRegister.h)

    set(application_files application/hoImageRegContainer2DRegistration.h)

    if (BUILD_CPU_OPTIMAL_FLOW_REG)

        add_library(pingvin_toolbox_cpureg SHARED
                ${opticalflow_files}
                ${transformation_files}
                ${solver_files}
                ${warper_files}
                ${similarity_files}
                ${register_files}
                ${application_files})

        set_target_properties(pingvin_toolbox_cpureg PROPERTIES VERSION ${PINGVIN_VERSION_STRING} SOVERSION ${PINGVIN_SOVERSION})


        target_link_libraries(pingvin_toolbox_cpureg
                pingvin_toolbox_operator
                pingvin_toolbox_mri_image
                pingvin_toolbox_cpucore
                pingvin_toolbox_cpu_solver
                pingvin_toolbox_log
                pingvin_toolbox_cpucore_math
                pingvin_toolbox_image_analyze_io
                armadillo
                )
        target_include_directories(pingvin_toolbox_cpureg
                PUBLIC

                $<INSTALL_INTERFACE:${PINGVIN_INSTALL_INCLUDE_PATH}>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/transformation>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/application>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dissimilarity>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/register>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solver>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/warper>
                )

    else ()
        message("cpureg library will not be built ... ")
    endif ()

    if (BUILD_CPU_OPTIMAL_FLOW_REG)
        source_group(opticalflow FILES ${opticalflow_files})
        install(FILES ${opticalflow_files} DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)
    endif ()

    source_group(transformation FILES ${transformation_files})
    source_group(solver FILES ${solver_files})
    source_group(warper FILES ${warper_files})
    source_group(similarity FILES ${similarity_files})
    source_group(register FILES ${register_files})
    source_group(application FILES ${application_files})

    if (BUILD_CPU_OPTIMAL_FLOW_REG)

		install(TARGETS pingvin_toolbox_cpureg
            EXPORT pingvin-export
			LIBRARY DESTINATION lib
			ARCHIVE DESTINATION lib
			RUNTIME DESTINATION bin
			COMPONENT main
			)

    endif ()

    install(FILES
            ${transformation_files}
            ${solver_files}
            ${warper_files}
            ${similarity_files}
            ${register_files}
            ${application_files}
            DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)