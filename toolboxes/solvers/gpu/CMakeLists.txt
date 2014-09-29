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
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  )

cuda_add_library(gadgetron_toolbox_gpusolvers SHARED 
    gpusolvers_export.h
    cuSolverUtils.cu
  )

target_link_libraries(gadgetron_toolbox_gpusolvers 
  gadgetron_toolbox_gpucore 
  ${Boost_LIBRARIES}
  ${CUDA_LIBRARIES}
  ${CUDA_CUBLAS_LIBRARIES} 
  )

install(TARGETS gadgetron_toolbox_gpusolvers DESTINATION lib)

install(FILES 	
  cuSbCgSolver.h
  cuSbcCgSolver.h
  cuCgPreconditioner.h
  cuLwSolver.h
  cuLbfgsSolver.h
  cuSbLwSolver.h
  cuSbcLwSolver.h
  cuCgSolver.h
  cuNlcgSolver.h
  cuGpBbSolver.h
  hoCuCgSolver.h
  hoCuNlcgSolver.h
  hoCuSbcCgSolver.h
  hoCuGpBbSolver.h
  cuSolverUtils.h
  gpusolvers_export.h
  DESTINATION include)
