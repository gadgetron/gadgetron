include(GoogleTest)

    if (WIN32)
        link_directories(${Boost_LIBRARY_DIRS})
    endif ()
    set(test_src_files
            tests.cpp
            hoNDArray_elemwise_test.cpp
            hoNDArray_blas_test.cpp
            hoNDArray_utils_test.cpp
            hoNDArray_reductions_test.cpp
            read_writer_test.cpp
            hoNDFFT_test.cpp
            hoNFFT_test.cpp
            hoNDWavelet_test.cpp
            curveFitting_test.cpp
            image_morphology_test.cpp
            IsmrmrdContextVariables_test.cpp
            StorageSpaces_test.cpp
            pattern_recognition_test.cpp
            cmr_mapping_test.cpp
            hoNDArray_linalg_test.cpp
            core_test.cpp
            core_primitive_io_test.cpp 
            threadpool_test.cpp
            from_string_test.cpp
            hoNDArrayView_test.cpp
            ChannelAlgorithmsTest.cpp
            cmr_strain_test.cpp
            cmr_thickening_test.cpp
            cmr_analytical_strain_test.cpp
            #lapack_test.cpp
            hoSDC_test.cpp
            nhlbi_compression_tests.cpp
            gadgets/setup_gadget.h 
            gadgets/AcquisitionAccumulateTrigget_test.cpp 
            gadgets/FlagTriggerParsing_test.cpp  
            )

    if (PYTHONLIBS_FOUND)
        set(test_src_files ${test_src_files} python_converter_test.cpp)
    endif ()

    if (CUDA_FOUND)
        set(test_src_files ${test_src_files}
                cuNDArray_elemwise_test.cpp
                cuNDArray_operators_test.cpp
                cuNDArray_blas_test.cpp
                cuNDArray_utils_test.cpp
                #        vector_td_test.cpp
                cuVector_td_test_kernels.h
                cuVector_td_test_kernels.cu
                cuNDFFT_test.cpp
                cuSDC_test.cpp
                )

        cuda_add_executable(test_all
                ${test_src_files}
                )
    else ()
        add_executable(test_all
                ${test_src_files})
    endif ()

    if (CUDA_FOUND)
        target_link_libraries(test_all
                gadgetron_toolbox_gpu
                )
    endif ()
    target_link_libraries(test_all
            gadgetron_core
            gadgetron_core_readers
            gadgetron_core_writers
            gadgetron_mricore
            gadgetron_toolbox_cpucore
            gadgetron_toolbox_cpucore_math
            gadgetron_toolbox_cpufft
            gadgetron_toolbox_cpunfft
            gadgetron_toolbox_cpudwt
            gadgetron_toolbox_cpu_image
            gadgetron_toolbox_cpucore
            gadgetron_toolbox_log
            gadgetron_toolbox_cpuklt
            gadgetron_toolbox_image_analyze_io
            gadgetron_toolbox_mri_core
            gadgetron_toolbox_cpuoperator
            gadgetron_toolbox_cpu_image
            gadgetron_toolbox_cmr
            gadgetron_toolbox_pr
            gadgetron_toolbox_cpusdc
            ${GTEST_LIBRARIES}
            gmock
            )

    if (PYTHONLIBS_FOUND)
        target_link_libraries(test_all
                gadgetron_toolbox_python
                python)
    endif ()

    gtest_discover_tests(test_all)

    install(TARGETS test_all DESTINATION bin COMPONENT main)

