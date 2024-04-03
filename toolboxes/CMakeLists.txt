if (WIN32)
    add_definitions(-DTIXML_USE_STL)
endif ()


add_subdirectory(log)

add_subdirectory(operators)
add_subdirectory(solvers)

add_subdirectory(fft)
add_subdirectory(core)
add_subdirectory(mri_core)
add_subdirectory(mri)
add_subdirectory(cmr)

add_subdirectory(mri_image)

add_subdirectory(klt)
add_subdirectory(fatwater)
add_subdirectory(deblur)
add_subdirectory(nfft)
add_subdirectory(dwt)
add_subdirectory(registration)

add_subdirectory(ffd)
add_subdirectory(image)
add_subdirectory(pattern_recognition)
add_subdirectory(denoise)

add_subdirectory(image_io)
add_subdirectory(T1)
# Should we compile the python toolbox
if (BUILD_PYTHON_SUPPORT)
    message("Find python; build gadgetron toolbox python")
    add_subdirectory(python)
else ()
    message("Cannot find python; do not build gadgetron toolbox python")
endif ()

if (ARMADILLO_FOUND AND Boost_FOUND AND PLPLOT_FOUND)
    add_subdirectory(plplot)
endif ()

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)

if(CUDA_FOUND)
    # We compile all GPU code into one library to workaround Thrust and Cub using statics in templates
    # (e.g. https://github.com/NVIDIA/cub/blob/48768e86a40c25d231a6ff1e04107a60016314a7/cub/util_device.cuh#L208)
    # that when used from multiple libraries can cause crashes.
    # See https://github.com/NVIDIA/thrust/issues/1401 for more information.
    if (WIN32)
        add_definitions(-D__BUILD_GADGETRON_GPU__)
        add_definitions(-D_USE_MATH_DEFINES)
        add_definitions(-DWIN32_LEAN_AND_MEAN)

        link_directories(${Boost_LIBRARY_DIRS})
    endif ()

    cuda_add_library(gadgetron_toolbox_gpu SHARED 
        core/gpu/check_CUDA.h
        core/gpu/CUBLASContextProvider.cpp
        core/gpu/CUBLASContextProvider.h
        core/gpu/cudaDeviceManager.cpp
        core/gpu/cudaDeviceManager.h
        core/gpu/cuNDArray_blas.cu
        core/gpu/cuNDArray_blas.h
        core/gpu/cuNDArray_elemwise.cu
        core/gpu/cuNDArray_elemwise.h
        core/gpu/cuNDArray_fileio.h
        core/gpu/cuNDArray_operators.cu
        core/gpu/cuNDArray_operators.h
        core/gpu/cuNDArray_reductions.cu
        core/gpu/cuNDArray_reductions.h
        core/gpu/cuNDArray_utils.cu
        core/gpu/cuNDArray_utils.h
        core/gpu/cuNDArray.h
        core/gpu/cuSparseMatrix.cu
        core/gpu/GadgetronCuException.h
        core/gpu/GPUTimer.h
        core/gpu/hoCuNDArray_blas.cpp
        core/gpu/hoCuNDArray_blas.h
        core/gpu/hoCuNDArray_elemwise.h
        core/gpu/hoCuNDArray_utils.h
        core/gpu/hoCuNDArray.h
        core/gpu/radial_utilities.cu
        core/gpu/radial_utilities.h
        core/gpu/real_utilities_device.h
        core/gpu/setup_grid.h
        dwt/gpu/cuNDDWT.cu
        dwt/gpu/cuNDDWT.h
        fft/gpu/cuFFTCachedPlan.h
        fft/gpu/cuFFTPlan.h
        fft/gpu/cuFFTPlan.hpp
        fft/gpu/cuNDFFT.cpp
        fft/gpu/cuNDFFT.cu
        fft/gpu/cuNDFFT.h
        mri/hyper/CSI_utils.cu
        mri/hyper/CSI_utils.h
        mri/hyper/CSIOperator.cpp
        mri/hyper/CSIOperator.h
        mri/hyper/gadgetron_toolbox_hyper_export.h
        mri/pmri/gpu/b1_map_NIH_Souheil.cu
        mri/pmri/gpu/b1_map.cu
        mri/pmri/gpu/b1_map.h
        mri/pmri/gpu/cuBuffer.cpp
        mri/pmri/gpu/cuBuffer.h
        mri/pmri/gpu/cuCartesianSenseOperator.cu
        mri/pmri/gpu/cuCartesianSenseOperator.h
        mri/pmri/gpu/cuNonCartesianKtSenseOperator.cu
        mri/pmri/gpu/cuNonCartesianKtSenseOperator.h
        mri/pmri/gpu/cuNonCartesianSenseOperator.cu
        mri/pmri/gpu/cuNonCartesianSenseOperator.h
        mri/pmri/gpu/cuSenseBuffer.cpp
        mri/pmri/gpu/cuSenseBuffer.h
        mri/pmri/gpu/cuSenseBufferCg.cpp
        mri/pmri/gpu/cuSenseBufferCg.h
        mri/pmri/gpu/cuSenseOperator.cu
        mri/pmri/gpu/cuSenseOperator.h
        mri/pmri/gpu/cuSpiritBuffer.cpp
        mri/pmri/gpu/cuSpiritOperator.h
        mri/pmri/gpu/gpupmri_export.h
        mri/pmri/gpu/htgrappa.cpp
        mri/pmri/gpu/htgrappa.cu
        mri/pmri/gpu/htgrappa.h
        mri/pmri/gpu/sense_utilities.cu
        mri/pmri/gpu/sense_utilities.h
        mri/pmri/gpu/senseOperator.h
        mri/pmri/gpu/spirit_calibration.cu
        mri/pmri/gpu/spirit_calibration.h
        mri/pmri/gpu/trajectory_utils.cu
        mri/pmri/gpu/trajectory_utils.h
        mri/sdc/gpu/cuSDC.cu
        nfft/gpu/cuGriddingConvolution.cu
        nfft/gpu/cuGriddingConvolution.h
        nfft/gpu/cuNFFT.cu
        nfft/gpu/cuNFFT.h
        nfft/gpu/cuNFFTOperator.cpp
        operators/gpu/../generalOperator.h
        operators/gpu/../linearOperator.h
        operators/gpu/cuConvolutionOperator.cu
        operators/gpu/cuConvolutionOperator.h
        operators/gpu/cuLaplaceOperator.cu
        operators/gpu/cuLaplaceOperator.h
        operators/gpu/cuPartialDerivativeOperator.cu
        operators/gpu/cuPartialDerivativeOperator.h
        operators/gpu/cuPartialDerivativeOperator2.cu
        operators/gpu/cuTv1dOperator.cu
        operators/gpu/cuTv1dOperator.h
        operators/gpu/cuTvOperator.cu
        operators/gpu/cuTvOperator.h
        registration/optical_flow/gpu/cuCKOpticalFlowSolver.cu
        registration/optical_flow/gpu/cuHSOpticalFlowSolver.cu
        registration/optical_flow/gpu/cuLinearResampleOperator.cu
        registration/optical_flow/gpu/cuOpticalFlowSolver.cu
        registration/optical_flow/gpu/cuResampleOperator.cu
        solvers/gpu/cuSolverUtils.cu
        solvers/gpu/gpusolvers_export.h
        dwt/gpu/cuHaarWaveletOperator.h
        dwt/gpu/cuHaarWaveletOperator.cu
    )

    target_include_directories(gadgetron_toolbox_gpu
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/core/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/dwt/gpu/
        ${CMAKE_CURRENT_SOURCE_DIR}/fft/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/mri/hyper
        ${CMAKE_CURRENT_SOURCE_DIR}/mri/pmri/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/nfft/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/operators/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/registration/optical_flow/
        ${CMAKE_CURRENT_SOURCE_DIR}/registration/optical_flow/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/mri/sdc/gpu
        ${CMAKE_CURRENT_SOURCE_DIR}/solvers
        ${CMAKE_CURRENT_SOURCE_DIR}/solvers/gpu
    )

    set_target_properties(gadgetron_toolbox_gpu PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

    target_link_libraries(gadgetron_toolbox_gpu 
        gadgetron_toolbox_cpucore
        gadgetron_toolbox_cpucore_math
        gadgetron_toolbox_log
        gadgetron_toolbox_nfft
        gadgetron_toolbox_operator
        gadgetron_toolbox_sdc
        ${CUDA_LIBRARIES}
        ${CUDA_CUFFT_LIBRARIES} 
        ${CUDA_CUBLAS_LIBRARIES} 
        ${CUDA_CUSPARSE_LIBRARIES}
        ${MKL_LIBRARIES}
    )


    install(TARGETS gadgetron_toolbox_gpu
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        COMPONENT main
    )

    install(FILES
        core/gpu/check_CUDA.h
        core/gpu/CUBLASContextProvider.h
        core/gpu/cudaDeviceManager.h
        core/gpu/cuNDArray_blas.h
        core/gpu/cuNDArray_converter.h
        core/gpu/cuNDArray_elemwise.h
        core/gpu/cuNDArray_fileio.h
        core/gpu/cuNDArray_math.h
        core/gpu/cuNDArray_operators.h
        core/gpu/cuNDArray_reductions.h
        core/gpu/cuNDArray_utils.h
        core/gpu/cuNDArray.h
        core/gpu/cuSparseMatrix.h
        core/gpu/GadgetronCuException.h
        core/gpu/GPUTimer.h				
        core/gpu/hoCuNDArray_blas.h
        core/gpu/hoCuNDArray_elemwise.h
        core/gpu/hoCuNDArray_math.h
        core/gpu/hoCuNDArray_utils.h
        core/gpu/hoCuNDArray.h
        core/gpu/radial_utilities.h
        core/gpu/real_utilities_device.h
        core/gpu/setup_grid.h
        dwt/gpu/cuDWTOperator.h
        dwt/gpu/cuNDDWT.h
        fft/gpu/cuFFTPlan.h
        fft/gpu/cuFFTPlan.hpp
        fft/gpu/cuNDFFT.h
        fft/gpu/gpufft_export.h
        mri/hyper/CSI_utils.h
        mri/hyper/CSIOperator.h
        mri/pmri/gpu/b1_map.h
        mri/pmri/gpu/cuBuffer.h
        mri/pmri/gpu/cuCartesianSenseOperator.h
        mri/pmri/gpu/cuNonCartesianKtSenseOperator.h
        mri/pmri/gpu/cuNonCartesianSenseOperator.h
        mri/pmri/gpu/cuSenseBuffer.h
        mri/pmri/gpu/cuSenseBufferCg.h
        mri/pmri/gpu/cuSenseOperator.h
        mri/pmri/gpu/cuSpiritBuffer.h
        mri/pmri/gpu/cuSpiritOperator.h
        mri/pmri/gpu/gpupmri_export.h
        mri/pmri/gpu/htgrappa.h
        mri/pmri/gpu/sense_utilities.h
        mri/pmri/gpu/senseOperator.h
        mri/sdc/gpu/cuSDC.h
        nfft/gpu/cuGriddingConvolution.h 
        nfft/gpu/cuNFFT.h
        operators/gpu/cuConvolutionOperator.h
        operators/gpu/cuDiagonalOperator.h
        operators/gpu/cuDiagonalSumOperator.h
        operators/gpu/cuDownsampleOperator.h
        operators/gpu/cuFFTOperator.h
        operators/gpu/cuImageOperator.h
        operators/gpu/cuLaplaceOperator.h
        operators/gpu/cuPartialDerivativeOperator.h
        operators/gpu/cuPartialDerivativeOperator2.h
        operators/gpu/cuTv1dOperator.h
        operators/gpu/cuTvOperator.h
        operators/gpu/cuTvPicsOperator.h
        operators/gpu/cuUpsampleOperator.h
        operators/gpu/gpuoperators_export.h
        operators/gpu/hoCuEncodingOperatorContainer.h
        operators/gpu/hoCuIdentityOperator.h
        operators/gpu/hoCuOperator.h
        operators/gpu/hoCuPartialDerivativeOperator.h
        operators/gpu/hoCuTvOperator.h
        operators/gpu/hoCuTvPicsOperator.h
        registration/optical_flow/gpu/cuCGHSOFSolver.h
        registration/optical_flow/gpu/cuCKOpticalFlowSolver.h
        registration/optical_flow/gpu/cuHSOpticalFlowSolver.h
        registration/optical_flow/gpu/cuLinearResampleOperator.h
        registration/optical_flow/gpu/cuOpticalFlowSolver.h
        registration/optical_flow/gpu/cuResampleOperator.h
        registration/optical_flow/gpu/gpureg_export.h
        solvers/gpu/cuCgPreconditioner.h
        solvers/gpu/cuCgSolver.h
        solvers/gpu/cuGpBbSolver.h
        solvers/gpu/cuLbfgsSolver.h
        solvers/gpu/cuLwSolver.h
        solvers/gpu/cuNlcgSolver.h
        solvers/gpu/cuSbcCgSolver.h
        solvers/gpu/cuSbCgSolver.h
        solvers/gpu/cuSbcLwSolver.h
        solvers/gpu/cuSbLwSolver.h
        solvers/gpu/cuSolverUtils.h
        solvers/gpu/gpusolvers_export.h
        solvers/gpu/hoCuCgSolver.h
        solvers/gpu/hoCuGpBbSolver.h
        solvers/gpu/hoCuNlcgSolver.h
        solvers/gpu/hoCuSbcCgSolver.h
        dwt/gpu/cuHaarWaveletOperator.h
        DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

else()
    message("CUDA not found, not compiling gpu toolbox")
endif()