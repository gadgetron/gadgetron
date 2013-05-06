include_directories( 
  ${CMAKE_SOURCE_DIR}/toolboxes/core 
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
  ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
  ${CMAKE_SOURCE_DIR}/toolboxes/operators
  ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers
  ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
  ${CMAKE_SOURCE_DIR}/apps/matlab
  ${Boost_INCLUDE_DIR} 
  ${ARMADILLO_INCLUDE_DIRS}
)

#add_subdirectory(MRI)
add_subdirectory(denoising)
#add_subdirectory(deblurring)
add_subdirectory(registration)
