if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUREG__)
endif (WIN32)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CUDA_INCLUDE_DIRS}
)

cuda_add_library(gadgetron_toolbox_gpureg SHARED 
  cuOpticalFlowSolver.cu 
  cuHSOpticalFlowSolver.cu 
  cuCKOpticalFlowSolver.cu 
  cuResampleOperator.cu 
  cuLinearResampleOperator.cu
#  cuRegistration_utils.cu
  )

set_target_properties(gadgetron_toolbox_gpureg PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpureg 
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_log
  ${CUDA_CUDART_LIBRARY} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
  )

install(TARGETS gadgetron_toolbox_gpureg DESTINATION lib COMPONENT main)

install(FILES
  cuOpticalFlowSolver.h
  cuHSOpticalFlowSolver.h
  cuCKOpticalFlowSolver.h
  gpureg_export.h
  cuResampleOperator.h
  cuLinearResampleOperator.h
#  cuRegistration_utils.h
  cuCGHSOFSolver.h
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)
