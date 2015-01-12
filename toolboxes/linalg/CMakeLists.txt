find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(FFTW3 COMPONENTS single double threads)

message("FFTW3_LIBRARIES: ${FFTW3_LIBRARIES}") 
if (WIN32)
	ADD_DEFINITIONS(-D__BUILD_GADGETRON_LINALG__)
endif (WIN32)


include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core/cpu)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/core)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/linalg)
include_directories(${CMAKE_SOURCE_DIR}/toolboxes/fft/cpu)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${ISMRMRD_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_library(linalg SHARED linalg_export.h 
                          matrix_vector_op.h 
                          matrix_vector_op.cpp 
                          matrix_decomposition.h 
                          matrix_decomposition.cpp)

add_executable(linalg_test linalg_test.cpp)

target_link_libraries(linalg ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${FFTW3_LIBRARIES})
target_link_libraries(linalg_test ${FFTW3_LIBRARIES} linalg)
target_link_libraries(linalg_test gadgetron_toolbox_cpucore_math gadgetron_toolbox_cpufft gadgetron_toolbox_log)

install(TARGETS linalg DESTINATION lib COMPONENT main)
install(TARGETS linalg_test DESTINATION bin COMPONENT main)

install(FILES matrix_vector_op.h matrix_decomposition.h linalg_export.h DESTINATION include COMPONENT main)
