if (WIN32)
  ADD_DEFINITIONS(-D__BUILD_GADGETRON_GPUSOLVERS__)
endif (WIN32)

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CUDA_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/gpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  )

cuda_add_library(gpusolvers SHARED 
  cuGpBbSolver.cu
  )

target_link_libraries(gpusolvers 
  gpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gpusolvers DESTINATION lib)

install(FILES 	
  cuSbCgSolver.h
  cuSbcCgSolver.h
  cuCgPreconditioner.h
  cuLwSolver.h
  cuSbLwSolver.h
  cuSbcLwSolver.h
  cuCgSolver.h
  cuGpBbSolver.h
  hoCuGpBbSolver.h
  gpusolvers_export.h
  DESTINATION include)
