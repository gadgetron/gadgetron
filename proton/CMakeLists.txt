find_package(HDF5 REQUIRED HL)


include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${GADGETRON_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/solvers
  ${CMAKE_SOURCE_DIR}/operators
)

cuda_add_library( proton_kernels SHARED
  cuOperatorPathBackprojection.cu hoCuOperatorPathBackprojection.cu proton_kernels.cu histogram.cu
)

target_link_libraries(proton_kernels gpucore hostutils ${CUDA_LIBRARIES})

add_executable(protonReconstruct protonReconstruct.cpp )

add_executable(hostProtonReconstruct hostProtonReconstruct.cpp)

add_executable(hostSubsetProtonReconstruct hostSubsetProtonReconstruct.cpp)

#add_executable(protonForward protonForwardProject.cpp)

#add_executable(LSQRTest LSQRTest.cpp)

target_link_libraries(protonReconstruct proton_kernels gpucore cpucore_math gpusolvers hostutils ${HDF5_LIBRARIES} ${CUDA_LIBRARIES} ${CULA_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(hostProtonReconstruct proton_kernels cpucore cpucore_math gpucore hostutils  ${HDF5_LIBRARIES}  ${CULA_LIBRARIES} ${CUDA_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(hostSubsetProtonReconstruct proton_kernels cpucore cpucore_math gpucore hostutils  ${HDF5_LIBRARIES}  ${CULA_LIBRARIES} ${CUDA_LIBRARIES} ${Boost_LIBRARIES})
#target_link_libraries(LSQRTest proton_kernels gpucore hostutils solvers operators ${CUDA_LIBRARIES})

#target_link_libraries(protonForward proton_kernels gpucore hostutils solvers operators ${CUDA_LIBRARIES})

install(TARGETS protonReconstruct DESTINATION bin)
#install(TARGETS hostProtonReconstruct DESTINATION bin)

#install(TARGETS protonForward DESTINATION bin)

