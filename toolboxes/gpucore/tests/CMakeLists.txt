
include_directories(${Boost_INCLUDE_DIR}
				${CMAKE_SOURCE_DIR}/toolboxes/gpucore
				${CMAKE_SOURCE_DIR}/toolboxes/cpucore)
				
add_library(gpucore_tests SHARED cuGTBLAS_test.cpp)
target_link_libraries(gpucore_tests gpucore)