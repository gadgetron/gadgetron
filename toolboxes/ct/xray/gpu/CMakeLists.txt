find_package(HDF5 REQUIRED HL)

if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUXRAY__)
endif (WIN32)


include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
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

cuda_add_library(gadgetron_toolbox_gpuxray SHARED
  conebeam_projection.cu 
  hoCuConebeamProjectionOperator.cpp 
  )

set_target_properties(gadgetron_toolbox_gpuxray PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(
  gadgetron_toolbox_gpuxray gadgetron_toolbox_gpucore gadgetron_toolbox_gpunfft
  gadgetron_toolbox_cpucore_math
  ${CUDA_LIBRARIES}
)

install(TARGETS gadgetron_toolbox_gpuxray DESTINATION lib COMPONENT main)

install(FILES 
  CBCT_acquisition.h
  CBCT_binning.h
  conebeam_projection.h
  hoCuConebeamProjectionOperator.h
  gpuxray_export.h 
  DESTINATION include COMPONENT main)
