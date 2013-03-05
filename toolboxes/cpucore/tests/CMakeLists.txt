
include_directories(${Boost_INCLUDE_DIR}
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore)
				
add_library(cpucore_tests SHARED hoGTBLAS_test.cpp hoNDArray_test.cpp)
target_link_libraries(cpucore_tests cpucore ${FFTW3_LIBRARIES} optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} ${ARMADILLO_LIBRARIES})