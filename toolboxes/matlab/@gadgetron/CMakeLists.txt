include_directories(${MATLAB_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/toolboxes/mri_core 
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
${CMAKE_SOURCE_DIR}/toolboxes/core/ 
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/operators
${CMAKE_SOURCE_DIR}/toolboxes/solvers
${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/matlab/ 
${ISMRMRD_INCLUDE_DIR}
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/ )

find_library(MATLAB_MEX_LIBRARY libmex.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)
find_library(MATLAB_MX_LIBRARY libmx.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)
find_library(MATLAB_ENG_LIBRARY libeng.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)


add_library(gadgetronNFFT SHARED gadgetronNFFT.cpp)
message(STATUS "Matlab libraries: ${MATLAB_LIBRARIES}")
target_link_libraries(
	gadgetronNFFT
    gadgetron_toolbox_matlab
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_gpunfft
    ${MATLAB_MEX_LIBRARY}
    ${MATLAB_MX_LIBRARY}
    ${MATLAB_ENG_LIBRARY}
	${CUDA_LIBRARIES} 
  	${CUDA_CUFFT_LIBRARIES} 
  	${CUDA_CUBLAS_LIBRARIES} 
    
)


add_library(gadgetronSENSE SHARED gadgetronSENSE.cpp)
target_link_libraries(
	gadgetronSENSE
    gadgetron_toolbox_matlab
    gadgetron_toolbox_cpucore
    gadgetron_toolbox_cpucore_math
    gadgetron_toolbox_gpunfft
    gadgetron_toolbox_gpuparallelmri
    ${MATLAB_MEX_LIBRARY}
    ${MATLAB_MX_LIBRARY}
    ${MATLAB_ENG_LIBRARY}
	${CUDA_LIBRARIES} 
  	${CUDA_CUFFT_LIBRARIES} 
  	${CUDA_CUBLAS_LIBRARIES} 
)


set(CMAKE_CXX_FLAGS "-O3 -DMX_COMPAT_32   -D_GNU_SOURCE  -DMATLAB_MEX_FILE  -DNDEBUG  -ansi -fPIC -fexceptions  -fno-omit-frame-pointer -std=c++11  -pthread  -Wl,--no-undefined  -Wl,-rpath,${MATLAB_LIB_PATHS}:${CMAKE_INSTALL_PREFIX}/lib:${CUDA_TOOLKIT_ROOT_DIR}/lib64" )


if (WIN32)
    if ( HAS_64_BIT )
        set(MATLAB_SUFFIX ".mexw64")
    else ()
        set(MATLAB_SUFFIX ".mexw32")
    endif ()

    add_definitions(-DMATLAB_DLL_EXPORT_SYM)
endif ()

if (UNIX)
    if ( HAS_64_BIT )
        set(MATLAB_SUFFIX ".mexa64")
    else ()
        set(MATLAB_SUFFIX ".mexglx")
    endif ()
endif ()

if (APPLE)
    if ( HAS_64_BIT )
        set(MATLAB_SUFFIX ".mexmaci64")
    endif ()
endif()


message("Matlab CXX flags: ${CMAKE_CXX_FLAGS}")
set_target_properties(gadgetronNFFT PROPERTIES PREFIX "" )
set_target_properties(gadgetronNFFT PROPERTIES SUFFIX ${MATLAB_SUFFIX} LINKER_LANGUAGE CXX)

set_target_properties(gadgetronSENSE PROPERTIES PREFIX "" )
set_target_properties(gadgetronSENSE PROPERTIES SUFFIX ${MATLAB_SUFFIX} LINKER_LANGUAGE CXX)


set(MATLAB_OUTPUT_DIR ${CMAKE_INSTALL_PREFIX}/matlab CACHE PATH "" )

install(TARGETS gadgetronNFFT DESTINATION ${MATLAB_OUTPUT_DIR} COMPONENT main)
install(TARGETS gadgetronSENSE DESTINATION ${MATLAB_OUTPUT_DIR} COMPONENT main)
