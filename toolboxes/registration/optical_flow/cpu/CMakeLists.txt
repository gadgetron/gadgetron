if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
endif (WIN32)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${ARMADILLO_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}/cpp
  ${ACE_INCLUDE_DIR}
)

add_library(cpureg ${LIBTYPE} 
    cpureg_export.h
    ../registrationSolver.h
    ../multiresRegistrationSolver.h
    ../resampleOperator.h
    ../opticalFlowSolver.h
    hoCKOpticalFlowSolver.h
    hoHSOpticalFlowSolver.h
    hoLinearResampleOperator.h
    hoLinearResampleOperator_eigen.h
    hoOpticalFlowSolver.h
    hoRegistration_utils.h
    hoOpticalFlowSolver.cpp
    hoHSOpticalFlowSolver.cpp
    hoCKOpticalFlowSolver.cpp
    hoLinearResampleOperator.cpp
    hoRegistration_utils.cpp
    )

target_link_libraries(cpureg 
  cpucore
  cpucore_math 
  ${ARMADILLO_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
  )

install(TARGETS cpureg DESTINATION lib)

install(FILES
    hoCKOpticalFlowSolver.h
    hoHSOpticalFlowSolver.h
    hoLinearResampleOperator.h
    hoLinearResampleOperator_eigen.h
    hoOpticalFlowSolver.h
    hoRegistration_utils.h
  DESTINATION include)

install(FILES hoLinearResampleOperator.h DESTINATION include)
