if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUREG__)
endif (WIN32)

include_directories(   
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
)

add_library(cpureg SHARED 
  hoOpticalFlowSolver.cpp
  hoHSOpticalFlowSolver.cpp
  hoCKOpticalFlowSolver.cpp
  hoLinearResampleOperator.cpp
  hoRegistration_utils.cpp
  )

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
