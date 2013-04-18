if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
endif (WIN32)

find_path(EIGEN_INCLUDE_DIR Eigen/Core
  HINTS /usr/include/eigen3 /usr/local/eigen3 /opt/local/include/eigen3)

IF(EIGEN_INCLUDE_DIR)
  message("Found Eigen. Exposing Eigen version of registration resampler.")
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})
ENDIF(EIGEN_INCLUDE_DIR)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
)

IF(EIGEN_INCLUDE_DIR)
  ADD_DEFINITIONS(-DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator.cpp
    hoLinearResampleOperator_eigen.cpp
    hoRegistration_utils.cpp
    )
ELSE(EIGEN_INCLUDE_DIR)
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator.cpp
    hoRegistration_utils.cpp
    )
ENDIF(EIGEN_INCLUDE_DIR)

target_link_libraries(cpureg 
  cpucore 
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS cpureg DESTINATION lib)

install(FILES
  hoOpticalFlowSolver.h
  hoHSOpticalFlowSolver.h
  hoCKOpticalFlowSolver.h
  hoLinearResampleOperator.h
  hoRegistration_utils.h
  DESTINATION include)

IF(EIGEN_INCLUDE_DIR)
  install(FILES hoLinearResampleOperator_eigen.h DESTINATION include)
ENDIF(EIGEN_INCLUDE_DIR)
