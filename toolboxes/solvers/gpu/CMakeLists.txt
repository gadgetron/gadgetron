if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_GPUSOLVERS__)
endif ()

if(WIN32)
  link_directories(${Boost_LIBRARY_DIRS})
endif()


cuda_add_library(gadgetron_toolbox_gpusolvers SHARED 
    gpusolvers_export.h
    cuSolverUtils.cu
  )

set_target_properties(gadgetron_toolbox_gpusolvers PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_toolbox_gpusolvers 
  gadgetron_toolbox_gpucore
  gadgetron_toolbox_solvers
  gadgetron_toolbox_log)

target_include_directories(gadgetron_toolbox_gpusolvers
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS gadgetron_toolbox_gpusolvers DESTINATION lib COMPONENT main)

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
  DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

set(GADGETRON_BUILD_RPATH "${CMAKE_CURRENT_BINARY_DIR};${GADGETRON_BUILD_RPATH}" PARENT_SCOPE)