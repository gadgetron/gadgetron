
include_directories( ${CUDA_INCLUDE_DIRS}
		     ${Boost_INCLUDE_DIR} 
		     ${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
		     ${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
		     ${CMAKE_SOURCE_DIR}/toolboxes/hostutils
		     ${CMAKE_SOURCE_DIR}/toolboxes/cpucore
		     ${CMAKE_SOURCE_DIR}/toolboxes/operators
		     ${CMAKE_SOURCE_DIR}/toolboxes/solvers)

add_executable(nfft main_nfft.cpp)
add_executable(nffth main_nffth.cpp)
add_executable(nffth_cg main_cg.cpp)
add_executable(nffth_sb main_sb.cpp)

target_link_libraries(nfft gpucore operators gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth gpucore operators  gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_cg solvers operators  gpucore gpunfft hostutils ${CUDA_LIBRARIES})
target_link_libraries(nffth_sb solvers operators  gpucore gpunfft hostutils ${CUDA_LIBRARIES})

install(TARGETS nfft nffth DESTINATION bin)
