if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUREG__)
endif (WIN32)

ADD_DEFINITIONS(-DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)

include_directories(   
  ${CUDA_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${EIGEN_INCLUDE_DIR}
)

cuda_add_library(gpureg SHARED 
  cuOpticalFlowSolver.cu 
  cuHSOpticalFlowSolver.cu 
  cuCKOpticalFlowSolver.cu 
  cuResampleOperator.cu 
  cuLinearResampleOperator.cu
  hoLinearResampleOperator.cpp
  )

target_link_libraries(gpureg gpucore ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})

install(TARGETS gpureg DESTINATION lib)

install(FILES
 registrationSolver.h
 multiresRegistrationSolver.h
 cuOpticalFlowSolver.h
 cuHSOpticalFlowSolver.h
 cuCKOpticalFlowSolver.h
 gpureg_export.h
 resampleOperator.h
 cuResampleOperator.h
 cuResampleOperator_macros.h
 cuLinearResampleOperator.h
 hoResampleOperator.h
 hoLinearResampleOperator.h
DESTINATION include)
