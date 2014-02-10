if ( HAS_64_BIT )

    if (WIN32)
        ADD_DEFINITIONS(-D__BUILD_GADGETRON_PLUS__)
    endif (WIN32)

    if(WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif(WIN32)

    if (MKL_FOUND)
        MESSAGE("MKL Found for gtPlus ... ")
        list(APPEND EXTRA_MKL_LIBRARIES mkl_core)
        if ( USE_OPENMP )
            list(APPEND EXTRA_MKL_LIBRARIES mkl_intel_thread)
        endif ( USE_OPENMP )

        INCLUDE_DIRECTORIES( ${MKL_INCLUDE_DIR} )
        LINK_DIRECTORIES( ${MKL_LIB_DIR} ${MKL_COMPILER_LIB_DIR} )
    endif (MKL_FOUND)

    include_directories(
        ${ACE_INCLUDE_DIR} 
        ${Boost_INCLUDE_DIR}
        ${ISMRMRD_INCLUDE_DIR}
        ${ISMRMRD_XSD_INCLUDE_DIR}
        ${XSD_INCLUDE_DIR}
        ${FFTW3_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${MKL_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core
        ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
        ${CMAKE_SOURCE_DIR}/toolboxes/operators
        ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers
        ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
        ${CMAKE_SOURCE_DIR}/gadgets/core
        ${HDF5_INCLUDE_DIR}
        ${HDF5_INCLUDE_DIR}/cpp
        ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
        ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
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
    )

    set( util_files util/gtPlusIOBase.h
        util/gtPlusIOBase.cpp
        util/gtPlusIOAnalyze.h
        util/gtPlusIOAnalyze.cpp
        util/gtPlusMemoryManager.h
        util/gtPlusMemoryManager.cpp 
        util/hoNDArrayMemoryManaged.h 
        util/hoNDArrayMemoryManaged.hxx )

    set( workflow_files workflow/gtPlusISMRMRDReconWorkFlow.h
        workflow/gtPlusISMRMRDReconWorkFlowCartesian.h
        workflow/gtPlusISMRMRDReconWorkFlowCartesian2DT.h
        workflow/gtPlusISMRMRDReconWorkFlowCartesian3DT.h
        workflow/gtPlusISMRMRDReconUtil.h
        workflow/gtPlusISMRMRDReconUtil.cpp
        workflow/gtPlusISMRMRDReconUtil.hxx
        workflow/gtPlusISMRMRDReconWorkOrder.h
        workflow/gtPlusISMRMRDReconWorkOrder2DT.h
        workflow/gtPlusISMRMRDReconWorkOrder3DT.h
        workflow/gtPlusISMRMRDReconWorker.h
        workflow/gtPlusISMRMRDReconWorker2DT.h
        workflow/gtPlusISMRMRDReconWorker3DT.h
        workflow/gtPlusISMRMRDReconWorker2DTGRAPPA.h
        workflow/gtPlusISMRMRDReconWorker2DTSPIRIT.h
        workflow/gtPlusISMRMRDReconWorker2DTL1SPIRITNCG.h
        workflow/gtPlusISMRMRDReconWorker2DTNoAcceleration.h
        workflow/gtPlusISMRMRDReconWorker3DTGRAPPA.h
        workflow/gtPlusISMRMRDReconWorker3DTSPIRIT.h
        workflow/gtPlusISMRMRDReconWorker3DTL1SPIRITNCG.h
        workflow/gtPlusISMRMRDReconWorker3DTNoAcceleration.h
        workflow/gtPlusCloudScheduler.h
        workflow/gtPlusCloudScheduler.cpp )

    set( algorithm_files algorithm/gtPlusAlgorithmBase.h 
                        algorithm/gtPlusGRAPPA.h 
                        algorithm/gtPlusSPIRIT.h
                        algorithm/gtPlusOperator.h 
                        algorithm/gtPlusSPIRITOperator.h 
                        algorithm/gtPlusSPIRIT2DOperator.h 
                        algorithm/gtPlusSPIRIT3DOperator.h 
                        algorithm/gtPlusSPIRIT2DTOperator.h 
                        algorithm/gtPlusSPIRITNoNullSpaceOperator.h 
                        algorithm/gtPlusSPIRITNoNullSpace2DOperator.h 
                        algorithm/gtPlusSPIRITNoNullSpace3DOperator.h 
                        algorithm/gtPlusSPIRITNoNullSpace2DTOperator.h 
                        algorithm/gtPlusWaveletOperator.h 
                        algorithm/gtPlusWavelet2DOperator.h 
                        algorithm/gtPlusWavelet3DOperator.h 
                        algorithm/gtPlusWaveletNoNullSpace2DOperator.h 
                        algorithm/gtPlusWaveletNoNullSpace3DOperator.h 
                        algorithm/gtPlusDataFidelityOperator.h )

    set( solver_files solver/gtPlusSolver.h 
                        solver/gtPlusLinearSolver.h 
                        solver/gtPlusNonLinearSolver.h
                        solver/gtPlusLSQRSolver.h 
                        solver/gtPlusNCGSolver.h )

    set( gtCloud_files config/gtCloud/myCloud_2DT.txt 
                    config/gtCloud/myCloud_3DT.txt 
                    config/gtCloud/myCloud_2DT_DualLayer.txt
                    config/gtCloud/myCloud_2DT_DualLayer_FirstLayer.txt )

    set( config_files config/GadgetronProgram_gtPlus_2DT_Cartesian.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_L1SPIRIT.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_SPIRIT.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_CloudNode.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_DualLayer_Gateway_L1SPIRIT.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_DualLayer_Gateway_SPIRIT.xml
                    config/GadgetronProgram_gtPlus_2DT_Cartesian_FirstLayer_CloudNode.xml
                    config/GadgetronProgram_gtPlus_2DT_FatWater.xml
                    config/GadgetronProgram_gtPlus_2DT_LGE.xml
                    config/GadgetronProgram_gtPlus_2DT_MOLLI.xml
                    config/GadgetronProgram_gtPlus_2DT_Perfusion.xml
                    config/GadgetronProgram_gtPlus_2DT_RealTimeCine.xml
                    config/GadgetronProgram_gtPlus_2DT_RealTimeFlow.xml
                    config/GadgetronProgram_gtPlus_2DT_T2W.xml
                    config/GadgetronProgram_gtPlus_3DT_Cartesian.xml
                    config/GadgetronProgram_gtPlus_3DT_Cartesian_CloudNode.xml
                    config/GadgetronProgram_gtPlus_3DT_Cartesian_L1SPIRIT.xml
                    config/GadgetronProgram_gtPlus_3DT_Cartesian_SPIRIT.xml )

    set( schema_files )

    # matlab
    if (MATLAB_FOUND)
        message("MATLAB FOUND: ${MATLAB_INCLUDE_DIR}, Matlab gt interface is being compiled.")
        SET(CMAKE_DEBUG_POSTFIX)
        include_directories( ${MATLAB_INCLUDE_DIR} )
        set( matlab_files matlab/gtMatlabConverter.h
                          matlab/gtMatlabConverterComplex.h )

    else(MATLAB_FOUND)
        message("MATLAB NOT FOUND ...")
        set( matlab_files )
    endif(MATLAB_FOUND)

    add_library(gtplus ${LIBTYPE} GtPlusExport.h ${util_files} ${workflow_files} ${algorithm_files} ${solver_files} ${config_files} ${schema_files} ${matlab_files})

    source_group(util FILES ${util_files})
    source_group(workflow FILES ${workflow_files})

    set( config_gtCloud_files ${config_files} ${gtCloud_files} )
    source_group(config FILES ${config_gtCloud_files})

    source_group(schema FILES ${schema_files})
    source_group(algorithm FILES ${algorithm_files})
    source_group(solver FILES ${solver_files})

    if (MATLAB_FOUND)
        source_group(matlab FILES ${matlab_files})
    endif(MATLAB_FOUND)

    target_link_libraries(gtplus cpucore cpucore_math ${MKL_LIBRARIES} ${EXTRA_MKL_LIBRARIES})

    if (CUDA_FOUND)
        target_link_libraries(gtplus gpuparallelmri)
    endif (CUDA_FOUND)

    install(TARGETS gtplus DESTINATION lib)

    # install gtplus files
    install (FILES  GtPlusExport.h
                    util/gtPlusIOBase.h 
                    util/gtPlusIOAnalyze.h 
                    util/hoNDArrayMemoryManaged.h 
                    util/hoNDArrayMemoryManaged.hxx 
                    util/gtPlusMemoryManager.h 
                    workflow/gtPlusISMRMRDReconWorkFlow.h
                    workflow/gtPlusISMRMRDReconWorkFlowCartesian.h
                    workflow/gtPlusISMRMRDReconWorkFlowCartesian2DT.h
                    workflow/gtPlusISMRMRDReconWorkFlowCartesian3DT.h
                    workflow/gtPlusISMRMRDReconUtil.h
                    workflow/gtPlusISMRMRDReconUtil.hxx
                    workflow/gtPlusISMRMRDReconWorkOrder.h
                    workflow/gtPlusISMRMRDReconWorker.h
                    workflow/gtPlusISMRMRDReconWorker2DT.h
                    workflow/gtPlusISMRMRDReconWorker3DT.h
                    workflow/gtPlusISMRMRDReconWorker2DTGRAPPA.h
                    workflow/gtPlusISMRMRDReconWorker2DTNoAcceleration.h
                    workflow/gtPlusISMRMRDReconWorker3DTGRAPPA.h
                    workflow/gtPlusISMRMRDReconWorker3DTNoAcceleration.h
                    algorithm/gtPlusAlgorithmBase.h
                    algorithm/gtPlusGRAPPA.h
                    algorithm/gtPlusSPIRIT.h
                    DESTINATION include)

    # install gadgetron program files
    install (FILES  ${config_files} 
                    DESTINATION config)

    install (FILES  ${gtCloud_files} 
                    DESTINATION config/gtCloud)

    install (FILES  ${schema_files} 
                    DESTINATION schema)

    if (MATLAB_FOUND)
        add_subdirectory(matlab)
    endif(MATLAB_FOUND)

    if (GTEST_FOUND)
        add_subdirectory(ut)
    endif (GTEST_FOUND)

endif ( HAS_64_BIT )
