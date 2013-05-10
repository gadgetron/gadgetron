if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
endif (WIN32)

find_path(EIGEN_INCLUDE_DIR Eigen/Core
  HINTS /usr/include/eigen3 /usr/local/eigen3 /opt/local/include/eigen3)

# Eigen MUST be at least version 3.1.3.
# There is no easy way to determine if this is the case, 
# instead the user needs to explicitly define the cmake variable EIGEN_SPARSE_IS_STABLE if this is the case...

IF(EIGEN_INCLUDE_DIR)
  message("Found Eigen.") 
  IF(EIGEN_SPARSE_IS_STABLE)  
    message("Compiling Eigen versions of the resampling operators.")
    INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})
  ELSE(EIGEN_SPARSE_IS_STABLE)  
    message("EIGEN_SPARSE_IS_STABLE is not defined, hence not compiling the Eigen versions of the resampling operators.")
	option(EIGEN_SPARSE_IS_STABLE "Eigen is at least version 3.1.3" Off)
  ENDIF(EIGEN_SPARSE_IS_STABLE)  
ELSE(EIGEN_INCLUDE_DIR)
  message("Eigen not found, not compiling the Eigen versions of the resampling operators.") 
ENDIF(EIGEN_INCLUDE_DIR)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${ARMADILLO_INCLUDE_DIR}
)

IF(EIGEN_INCLUDE_DIR AND EIGEN_SPARSE_IS_STABLE)
IF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator.cpp
    hoLinearResampleOperator_eigen.cpp
    hoRegistration_utils.cpp
    )
ELSE(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator_eigen.cpp
    hoRegistration_utils.cpp
    )
ENDIF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
ELSE(EIGEN_INCLUDE_DIR AND EIGEN_SPARSE_IS_STABLE)
IF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator.cpp
    hoRegistration_utils.cpp
    )
ELSE(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  add_library(cpureg SHARED 
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoRegistration_utils.cpp
    )
ENDIF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
ENDIF(EIGEN_INCLUDE_DIR AND EIGEN_SPARSE_IS_STABLE)

target_link_libraries(cpureg 
  cpucore
  cpucore_math 
  ${ARMADILLO_LIBRARIES}
  )

install(TARGETS cpureg DESTINATION lib)

install(FILES
  hoOpticalFlowSolver.h
  hoHSOpticalFlowSolver.h
  hoCKOpticalFlowSolver.h
  hoRegistration_utils.h
  DESTINATION include)

IF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )
  install(FILES hoLinearResampleOperator.h DESTINATION include)
ENDIF(${ARMADILLO_VERSION_STRING} VERSION_GREATER "3.39" )

IF(EIGEN_INCLUDE_DIR AND EIGEN_SPARSE_IS_STABLE)
  install(FILES hoLinearResampleOperator_eigen.h DESTINATION include)
ENDIF(EIGEN_INCLUDE_DIR AND EIGEN_SPARSE_IS_STABLE)

