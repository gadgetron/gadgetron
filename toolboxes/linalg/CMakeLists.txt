if (WIN32)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_LINALG__)
endif (WIN32)

include_directories( ${Boost_INCLUDE_DIR}
                     ${CMAKE_SOURCE_DIR}/toolboxes/ndarray
                     ${CMAKE_SOURCE_DIR}/toolboxes/linalg )

add_library(linalg SHARED matrix_vector_op.cpp matrix_decomposition.cpp)
add_executable(linalg_test linalg_test.cpp)

target_link_libraries(linalg ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(linalg_test linalg)

install(TARGETS linalg DESTINATION lib)
install(TARGETS linalg_test DESTINATION bin)

install(FILES matrix_vector_op.h matrix_decomposition.h linalg_export.h DESTINATION include)
