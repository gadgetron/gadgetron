include_directories(
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
  ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
  )

find_path(EIGEN_INCLUDE_DIR Eigen/Core
  HINTS /usr/include/eigen3 /usr/local/eigen3 /opt/local/include/eigen3)

IF(EIGEN_INCLUDE_DIR)
  ADD_DEFINITIONS(-DEIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET)
  INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIR})
  add_subdirectory(2d)
  add_subdirectory(3d)
ENDIF(EIGEN_INCLUDE_DIR)

