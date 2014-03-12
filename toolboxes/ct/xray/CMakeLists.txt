find_package(HDF5 REQUIRED HL)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
)

cuda_add_library(gpuxray SHARED
  conebeam_projection.cu 
  hoCuConebeamProjectionOperator.cpp 
  )

target_link_libraries(
  gpuxray gpucore gpunfft
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
