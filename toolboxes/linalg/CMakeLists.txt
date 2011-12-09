if (WIN32)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_LINALG__)
endif (WIN32)

include_directories(${CMAKE_SOURCE_DIR}/toolboxes/ndarray)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/linalg)

add_library(linalg SHARED matrix_vector_op.cpp)
add_executable(linalg_test linalg_test.cpp)

target_link_libraries(linalg ${BLAS_LIBRARIES})
target_link_libraries(linalg_test linalg)

install(TARGETS linalg DESTINATION lib)
install(TARGETS linalg_test DESTINATION bin)

install(FILES matrix_vector_op.h linalg_export.h DESTINATION include)
