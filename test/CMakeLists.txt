ENABLE_TESTING()

if(WIN32)
link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories(
  ${CMAKE_SOURCE_DIR}/operators
  ${Boost_INCLUDE_DIR}
  ${CUDA_INCLUDE_DIRS}
  ${GTEST_INCLUDE_DIRS}
  )
cuda_add_executable(test_all 
  tests.cpp 
  waveletTest.cpp
  )

target_link_libraries(test_all 
  cpucore 
  cpucore_math
  gpucore
  tomography_operators
  ${BOOST_LIBRARIES}	
  ${GTEST_LIBRARIES} 
  ${CUDA_LIBRARIES} 
  )

add_test(test_all test_all)
