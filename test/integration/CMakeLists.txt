
# detect the environment
# CUDA_FOUND, MKL_FOUND, WIN32, UNIX, APPLE, PYTHONLIBS_FOUND

# read in what will be tested for any environment
FILE(GLOB TEST_CASE_ALWAYS RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/always/*.cfg  )
FILE(GLOB TEST_CASE_WITH_MKL RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/mkl/*.cfg  )
FILE(GLOB TEST_CASE_WITH_GPU RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/gpu/*.cfg  )
FILE(GLOB TEST_CASE_WITH_PYTHON RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/python/*.cfg  )

message("Test cases found : ${TEST_CASE_ALWAYS}")
message("Test cases with MKL found : ${TEST_CASE_WITH_MKL}")
message("Test cases with_GPU : ${TEST_CASE_WITH_GPU}")
message("Test cases with python : ${TEST_CASE_WITH_PYTHON}")

# according to the building environment, generate the list of available test cases

set(TEST_CASE ${TEST_CASE_ALWAYS})

if (MKL_FOUND)
    set(TEST_CASE ${TEST_CASE} ${TEST_CASE_WITH_MKL})
endif (MKL_FOUND)

if (CUDA_FOUND)
    set(TEST_CASE ${TEST_CASE} ${TEST_CASE_WITH_GPU})
endif (CUDA_FOUND)

if (PYTHONLIBS_FOUND)
    set(TEST_CASE ${TEST_CASE} ${TEST_CASE_WITH_PYTHON})
endif (PYTHONLIBS_FOUND)

message("Test cases to be performed : ${TEST_CASE}")

# write the available test cases to a file
file(WRITE ./test_cases.txt)
foreach(test ${TEST_CASE})
    file(APPEND ./test_cases.txt "${test}\n")
endforeach()