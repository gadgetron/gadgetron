include_directories(${GTEST_INCLUDE})

add_executable(test_all tests.cpp)
target_link_libraries(test_all cpucore cpucore_tests gpucore gpucore_tests ${GTEST_LIBRARIES} gpunfft hostutils ${CUDA_LIBRARIES})

add_test(test_all test_all)
