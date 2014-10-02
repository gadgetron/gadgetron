find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


if (WIN32)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_LINALG__)
endif (WIN32)


include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core/cpu)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/linalg)
include_directories(${Boost_INCLUDE_DIR})

add_library(linalg SHARED matrix_vector_op.cpp matrix_decomposition.cpp)
add_executable(linalg_test linalg_test.cpp)

target_link_libraries(linalg ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(linalg_test linalg gadgetron_toolbox_cpucore_math)

install(TARGETS linalg DESTINATION lib)
install(TARGETS linalg_test DESTINATION bin)

install(FILES matrix_vector_op.h matrix_decomposition.h linalg_export.h DESTINATION include)
