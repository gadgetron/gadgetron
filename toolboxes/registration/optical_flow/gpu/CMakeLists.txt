if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUREG__)
endif (WIN32)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CUDA_INCLUDE_DIRS}
)

cuda_add_library(gpureg SHARED 
  cuOpticalFlowSolver.cu 
  cuHSOpticalFlowSolver.cu 
  cuCKOpticalFlowSolver.cu 
  cuResampleOperator.cu 
  cuLinearResampleOperator.cu
  cuRegistration_utils.cu
  )

target_link_libraries(gpureg 
  gpucore 
  ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
  )

install(TARGETS gpureg DESTINATION lib)

install(FILES
  cuOpticalFlowSolver.h
  cuHSOpticalFlowSolver.h
  cuCKOpticalFlowSolver.h
  gpureg_export.h
  cuResampleOperator.h
  cuLinearResampleOperator.h
  cuRegistration_utils.h
  DESTINATION include)
