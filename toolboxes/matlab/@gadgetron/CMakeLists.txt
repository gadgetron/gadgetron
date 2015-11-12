include_directories(${MATLAB_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/toolboxes/mri_core 
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math 
${CMAKE_SOURCE_DIR}/toolboxes/core/ 
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/core/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/operators
${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu 
${CMAKE_SOURCE_DIR}/toolboxes/matlab/ 
${ISMRMRD_INCLUDE_DIR}
${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/ )

FIND_LIBRARY(MATLAB_MEX_LIBRARY libmex.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)
FIND_LIBRARY(MATLAB_MX_LIBRARY libmx.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)
FIND_LIBRARY(MATLAB_ENG_LIBRARY libeng.so ${MATLAB_LIB_PATHS} NO_DEFAULT_PATH)


add_library(gadgetronNFFT SHARED gadgetronNFFT.cpp)
MESSAGE(STATUS "Matlab libraries: ${MATLAB_LIBRARIES}")
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


set(CMAKE_CXX_FLAGS "-O3 -DMX_COMPAT_32   -D_GNU_SOURCE  -DMATLAB_MEX_FILE  -DNDEBUG  -ansi -fPIC -fexceptions  -fno-omit-frame-pointer -std=c++11  -pthread  -Wl,--no-undefined -Wl,-soname,gadgetronNFFT.mexa64 -Wl,-rpath,${MATLAB_LIB_PATHS}:${CMAKE_INSTALL_PREFIX}/lib:${CUDA_TOOLKIT_ROOT_DIR}/lib64" )

MESSAGE("Matlab CXX flags: ${CMAKE_CXX_FLAGS}")
set_target_properties(gadgetronNFFT PROPERTIES PREFIX "" )
set_target_properties(gadgetronNFFT PROPERTIES SUFFIX .mexa64 LINKER_LANGUAGE CXX)

SET(MATLAB_OUTPUT_DIR PATH CACHE PATH "")

install(DIRECTORY DESTINATION ${MATLAB_OUTPUT_DIR}/@gadgetron)
install(TARGETS gadgetronNFFT DESTINATION ${MATLAB_OUTPUT_DIR}/@gadgetron COMPONENT main)
