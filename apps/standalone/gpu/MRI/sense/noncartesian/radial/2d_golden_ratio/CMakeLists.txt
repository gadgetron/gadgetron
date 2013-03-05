include_directories(    ${CUDA_INCLUDE_DIRS}
				        ${Boost_INCLUDE_DIR} 
						${CMAKE_SOURCE_DIR}/toolboxes/gpucore 
						${CMAKE_SOURCE_DIR}/toolboxes/gpunfft
						${CMAKE_SOURCE_DIR}/toolboxes/hostutils
						${CMAKE_SOURCE_DIR}/toolboxes/gpuct
						${CMAKE_SOURCE_DIR}/toolboxes/solvers
						${CMAKE_SOURCE_DIR}/toolboxes/operators
						${CMAKE_SOURCE_DIR}/toolboxes/gpupmri
						${CMAKE_SOURCE_DIR}/toolboxes/cpucore)

add_executable(radial_sense_cg main_cg.cpp)
add_executable(radial_sense_sb main_sb.cpp)
add_executable(radial_sense_sbc main_sbc.cpp)

target_link_libraries(radial_sense_cg operators gpucore gpuparallelmri gpunfft hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(radial_sense_sb operators gpucore gpuparallelmri gpunfft hostutils solvers ${CUDA_LIBRARIES})
target_link_libraries(radial_sense_sbc operators gpucore gpuparallelmri gpunfft hostutils solvers ${CUDA_LIBRARIES})

install(TARGETS radial_sense_cg radial_sense_sb radial_sense_sbc DESTINATION bin)
