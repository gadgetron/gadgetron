find_package(HDF5 REQUIRED HL)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}/cpp
  ${ISMRMRD_INCLUDE_DIR}
  ${ARMADILLO_INCLUDE_DIRS}
)

cuda_add_library(gpuxray SHARED
  conebeam_projection.cu 
  hoCuConebeamProjectionOperator.cpp 
  )

target_link_libraries(
  gpuxray gpucore gpunfft
  cpucore_math
  ${CUDA_LIBRARIES}
)

install(TARGETS gpuxray DESTINATION lib)

install(FILES 
  CBCT_acquisition.h
  CBCT_binning.h
  conebeam_projection.h
  hoCuConebeamProjectionOperator.h
  gpuxray_export.h 
  DESTINATION include)
