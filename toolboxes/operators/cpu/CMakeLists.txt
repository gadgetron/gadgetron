if (WIN32)
ADD_DEFINITIONS(-D__BUILD_GADGETRON_CPUOPERATORS__)
endif (WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  )

#target_link_libraries(cpuoperators  
#  )

#install(TARGETS cpuoperators DESTINATION lib)

install(FILES 	
  hoIdentityOperator.h
  hoImageOperator.h
  hoDiagonalOperator.h
  hoDiagonalSumOperator.h
  hoFFTOperator.h
  hoPartialDerivativeOperator.h
  hoTvOperator.h
  hoTvPicsOperator.h
  DESTINATION include)
