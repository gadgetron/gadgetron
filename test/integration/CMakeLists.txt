
# detect the environment
# CUDA_FOUND, MKL_FOUND, WIN32, UNIX, APPLE, PYTHONLIBS_FOUND

# read in what will be tested for any environment
FILE(GLOB TEST_CASE_ALWAYS RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/*.cfg  )

set(TEST_CASE ${TEST_CASE_ALWAYS})

# Remove this after WIN python issue is resolved
if (WIN32)
  string(REPLACE "simple_gre_python" "simple_gre" TEST_CASE2 ${TEST_CASE})
  string(REPLACE "cfg" "cfg;" TEST_CASE ${TEST_CASE2})
endif (WIN32)

# write the available test cases to a file
file(WRITE ./test_cases.txt)
foreach(test ${TEST_CASE})
    file(APPEND ./test_cases.txt "${test}\n")
endforeach()
