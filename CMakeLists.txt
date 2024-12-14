cmake_minimum_required(VERSION 3.17...3.21 FATAL_ERROR)

if (${CMAKE_VERSION} VERSION_LESS 3.17)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else ()
    cmake_policy(VERSION 3.17)
endif()
project(PINGVIN LANGUAGES CXX C)

find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif (CCACHE_FOUND)

if (EXISTS $ENV{CMAKE_PREFIX_PATH})
    list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif ()
if (EXISTS $ENV{TORCH_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{TORCH_HOME})
endif ()
if (CMAKE_PREFIX_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
endif ()
message(STATUS "Looking for packages in : ${CMAKE_PREFIX_PATH}")
# ---

# check the compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (NOT CMAKE_CXX_COMPILER_VERSION) #Compiler version is not set on Ubuntu 12.02 (gcc 4.6)
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
    endif ()
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message(STATUS "CXX COMPILER VERSION: ${CMAKE_CXX_COMPILER_ID} : ${CMAKE_CXX_COMPILER_VERSION}")
        message(FATAL_ERROR "Pingvin requires GCC version >= 6.0")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        message(FATAL_ERROR "Pingvin requires Clang version >= 3.5")
    endif ()
    add_compile_definitions(BOOST_NO_CXX98_FUNCTION_BASE)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        message(FATAL_ERROR "Pingvin requires MSVC 2013 or above ")
    endif ()
else ()
    message(WARNING "Unsupported compiler!")
endif ()

#VERSIONING AND INSTALL PATHS
set(PINGVIN_VERSION_MAJOR 0)
set(PINGVIN_VERSION_MINOR 1)
set(PINGVIN_VERSION_PATCH 0)
set(PINGVIN_VERSION_STRING ${PINGVIN_VERSION_MAJOR}.${PINGVIN_VERSION_MINOR}.${PINGVIN_VERSION_PATCH})
set(PINGVIN_SOVERSION ${PINGVIN_VERSION_MAJOR}.${PINGVIN_VERSION_MINOR})
find_package(Git)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE PINGVIN_GIT_SHA1 ERROR_VARIABLE PINGVIN_GIT_STDERR)
    string(STRIP "${PINGVIN_GIT_SHA1}" PINGVIN_GIT_SHA1)
    string(LENGTH "${PINGVIN_GIT_SHA1}" PINGVIN_GIT_SHA1_LEN)
    if (${PINGVIN_GIT_SHA1_LEN} LESS 40)
        message(WARNING "Could not determine SHA-1 hash: ${PINGVIN_GIT_STDERR}")
        set(PINGVIN_GIT_SHA1 "NA")
    endif ()
else ()
    set(PINGVIN_GIT_SHA1 "NA")
endif ()
set(PINGVIN_INSTALL_CMAKE_PATH share/pingvin/cmake)
set(PINGVIN_INSTALL_CONFIG_PATH share/pingvin/config)
set(PINGVIN_INSTALL_PYTHON_MODULE_PATH share/pingvin/python)
set(PINGVIN_INSTALL_SCHEMA_PATH share/pingvin/schema)
set(PINGVIN_INSTALL_INCLUDE_PATH include/pingvin)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# create the Pingvin version header
configure_file("${CMAKE_SOURCE_DIR}/cmake/pingvin_sha1.h.in" ${PROJECT_BINARY_DIR}/include/pingvin_sha1.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/include)
install(FILES ${PROJECT_BINARY_DIR}/include/pingvin_sha1.h DESTINATION ${PINGVIN_INSTALL_INCLUDE_PATH} COMPONENT main)


#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

# build options for 64 bits system
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(STATUS "64bit system is found")
    set(HAS_64_BIT On CACHE BOOL "64bit build")
else ()
    message(STATUS "32bit system is found")
    set(HAS_64_BIT Off CACHE BOOL "64bit build")
endif ()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# build options for OpenMP support
find_package(OpenMP)
option(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if (USE_OPENMP)
        message(STATUS "OpenMP multithreading enabled")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSE_OMP)
    else ()
        message(STATUS "OpenMP multithreading is supported, but disabled")
    endif ()
else ()
    message(STATUS "OpenMP multithreading not supported")
endif ()

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(CheckCXXCompilerFlag)
if (APPLE)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -Wno-array-bounds -Wno-switch -Wno-dynamic-class-memaccess -Wno-undefined-var-template -Wno-parentheses")
    endif ()
endif ()

add_compile_options(-Werror)
check_cxx_compiler_flag(-Werror=return-type COMPILER_SUPPORTS_WRETURNTYPE)
if (COMPILER_SUPPORTS_WRETURNTYPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
endif ()
check_cxx_compiler_flag(-Werror=terminate COMPILER_SUPPORTS_WTERMINATE)
if (COMPILER_SUPPORTS_WTERMINATE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=terminate")
endif ()
check_cxx_compiler_flag(-Werror=non-virtual-dtor COMPILER_SUPPORT_WNON_VIRTUAL)
if (COMPILER_SUPPORT_WNON_VIRTUAL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=non-virtual-dtor")
endif ()

check_cxx_compiler_flag(-Werror=format-overflow COMPILER_SUPPORT_WOVERFLOW)
if (COMPILER_SUPPORT_WOVERFLOW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=format-overflow")
endif()

# whether to suppress compilation warnings
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif ()

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# We actually only use system and thread explicitly, but they require linking in date_time and chrono
find_package(Boost 1.80.0 COMPONENTS coroutine system date_time program_options filesystem timer REQUIRED )
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS") # To suppress messages in later versions of boost, in version 1.76.0 there are internal conflicts in boost
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS") # To suppress messages in later versions of boost, in version 1.76.0 there are internal conflicts in boost
if (NOT DEFINED Boost_VERSION_STRING)
    set(Boost_VERSION_STRING "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
endif()

option(USE_CUDA "Build with CUDA support" On)
if (USE_CUDA)
    find_package(CUDA_advanced)
endif()
if (CUDA_FOUND AND USE_CUDA)
    if (NOT CUDA_cublas_device_LIBRARY)
        list(REMOVE_ITEM CUDA_CUBLAS_LIBRARIES CUDA_cublas_device_LIBRARY-NOTFOUND)
    endif()

    add_definitions(-DUSE_CUDA)
    set( PINGVIN_CUDA_FOUND_BOOL 1 )
    include_directories( ${CUDA_INCLUDE_DIRS} )
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --std=c++17")
    set(CMAKE_CUDA_STANDARD 17)

    # Compile kernels for compute models 1.0 and 2.0 as default for Cuda 4.1
    # Support compute model 3.0 from Cuda 4.2 and up
    # Support compute model 3.5 from Cuda 5 and up

    option(PINGVIN_CUDA_ALL_COMPUTE_MODEL "Build CUDA components for all computing models" OFF)
    if (PINGVIN_CUDA_ALL_COMPUTE_MODEL)
        message(STATUS "Compiling CUDA components to support all compute model 3.5, 5.0, 5.2, 6.0, 6.1, 7.0")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_35,code=sm_35")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_50,code=sm_50")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_52,code=sm_52")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_60,code=sm_60")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_61,code=sm_61")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_70,code=sm_70")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_75,code=sm_75")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_80,code=sm_80")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-gencode arch=compute_86,code=sm_86")
    endif ()

    message(STATUS "Compiling with ${CUDA_NVCC_FLAGS}")
    message(STATUS "CUDA libraries are ${CUDA_LIBRARIES}, ${CUDA_CUFFT_LIBRARIES}, ${CUDA_CUBLAS_LIBRARIES}, ${CUDA_CUSPARSE_LIBRARIES}")
else ()
    message(STATUS "CUDA not found. CUDA components will not be compiled.")
    set( PINGVIN_CUDA_FOUND_BOOL 0 )
    set( CUDA_NVCC_FLAGS "NA")
endif ()

option(REQUIRE_SIGNED_CONFIG "Force Pingvin to accept only signed config files." Off)
if (REQUIRE_SIGNED_CONFIG)
    message(STATUS "Compile with signed config; GTBabylon is needed")
    add_definitions(-DUSE_GTBABYLON)
    find_package(GTBabylon REQUIRED)
endif()

find_package(Armadillo REQUIRED)
add_library(armadillo INTERFACE)
set_property(TARGET armadillo PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ARMADILLO_INCLUDE_DIR})
add_definitions(-DARMA_DONT_USE_WRAPPER -DARMA_USE_CXX11 -DARMA_64BIT_WORD)
install(TARGETS armadillo EXPORT pingvin-export)

add_subdirectory(libmrd)

find_package(FFTW3 REQUIRED COMPONENTS single double)
add_library(FFTW INTERFACE)
target_link_libraries(FFTW INTERFACE ${FFTW3_LIBRARIES})
set_property(TARGET FFTW PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FFTW3_INCLUDE_DIR})
install(TARGETS FFTW EXPORT pingvin-export)

option(USE_MKL Off)
if (USE_MKL)
    if (LINUX)
        list(APPEND CMAKE_PREFIX_PATH "/opt/intel/oneapi/mkl/latest")
    endif()
    find_package(MKL REQUIRED)
    add_library(GTBLAS INTERFACE)
    set_property(TARGET GTBLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MKL_INCLUDE_DIR})
    set_property(TARGET GTBLAS PROPERTY INTERFACE_LINK_LIBRARIES ${MKL_LIBRARIES})

    add_definitions(-DUSE_MKL -DARMA_BLAS_LONG_LONG -DARMA_USE_BLAS -DARMA_USE_LAPACK)
else ()
    find_package(CBLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    add_library(GTBLAS INTERFACE)
    set_property(TARGET GTBLAS PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CBLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR} ${LAPACKE_INCLUDE_DIR})
    set_property(TARGET GTBLAS PROPERTY INTERFACE_LINK_LIBRARIES ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES})

    find_package(PkgConfig)
    pkg_check_modules(PC_lapacke lapacke)

    string(COMPARE EQUAL ${BLA_VENDOR} "OpenBLAS" USING_OPENBLAS)
    if (USING_OPENBLAS)
        message(STATUS "OpenBLAS detected; setting max number of OpenBLAS pthreads to one.")
        add_definitions(-DFORCE_LIMIT_OPENBLAS_NUM_THREADS)
    endif()
endif ()

install(TARGETS GTBLAS EXPORT pingvin-export)

if(EXISTS $ENV{TORCH_HOME})
    find_package(Torch)
    if (TORCH_FOUND)
        add_definitions(-DUSE_TORCH)
        message(STATUS "Found Pytorch c++ library at ${TORCH_LIBRARIES}")
    else()
        message(STATUS "Did not find Pytorch c++ library")
    endif ()
else()
    message(STATUS "Did not find TORCH_HOME")
endif()

option(BUILD_PYTHON_SUPPORT "Build Python gadgets etc" On)
if (BUILD_PYTHON_SUPPORT)
    find_package(Python3 COMPONENTS Development NumPy REQUIRED)
    find_package(Boost COMPONENTS python3${Python3_VERSION_MINOR} REQUIRED)

    message(STATUS "Found NumPy:")
    message(STATUS "NumPy include: ${Python3_NumPy_INCLUDE_DIRS}")
    message(STATUS "NumPy libraries: ${Python3_NumPy_LIBRARIES}")
    message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")

    if (Python3_VERSION_MINOR LESS 7)
        message(FATAL_ERROR "Pingvin Python support requires Python 3.7 or newer. Found: ${Python3_VERSION}")
    endif()

    message(STATUS "Boost Version: ${Boost_VERSION_STRING} (with Python ${Python3_VERSION})")

    add_library(python INTERFACE)
    set_property(TARGET python PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
    target_link_libraries(python INTERFACE Python3::Python Python3::NumPy Boost::python3${Python3_VERSION_MINOR})

    install(TARGETS python EXPORT pingvin-export)
endif()

option(BUILD_TESTING "Enable test building" On)
if (BUILD_TESTING)
    enable_testing()
    find_package(GTest)
    if (NOT GTEST_FOUND)
        find_path(GTEST_SRC_DIR src/gtest.cc HINTS /usr/src/gtest)
        find_path(GTEST_INCLUDE_DIRS gtest.h HINTS /usr/include/gtest)
        if (GTEST_SRC_DIR AND GTEST_INCLUDE_DIRS)
            message(STATUS "GTest src package found. Compiling as part of Pingvin.")
            add_subdirectory(${GTEST_SRC_DIR} ${CMAKE_BINARY_DIR}/gtest )
            include_directories(${GTEST_INCLUDE_DIRS})
            set(GTEST_FOUND 1)
            set(GTEST_LIBRARIES GTest::gtest GTest::gtest_main)
        else()
            message(FATAL_ERROR "GTest required for building tests")
        endif ()

    endif ()
    MESSAGE(STATUS "GTEST FOUND: ${GTEST_FOUND}")
    add_subdirectory(test)
else()
    MESSAGE(STATUS "Testing not being built")
endif()

add_subdirectory(apps)
add_subdirectory(cmake)
add_subdirectory(core)
add_subdirectory(gadgets)
add_subdirectory(toolboxes)

option(BUILD_WITH_GPERFTOOLS_PROFILER "Build with the gperftools profiler." Off)
if (BUILD_WITH_GPERFTOOLS_PROFILER)
    message("Adding gpertools to build.")
    find_package(Gperftools REQUIRED)
endif()

option(BUILD_DOCUMENTATION "Build Pingvin documentation." Off)
if(BUILD_DOCUMENTATION)
    add_subdirectory(doc)
endif()


#  ---   Main Library  (end) ----
include(GNUInstallDirs)
install(EXPORT pingvin-export
        FILE pingvin-targets.cmake
        NAMESPACE
        Pingvin::
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/pingvin)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/pingvin-config.cmake.in
    ${CMAKE_BINARY_DIR}/cmake/pingvin-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pingvin
    PATH_VARS PINGVIN_INSTALL_PYTHON_MODULE_PATH PINGVIN_INSTALL_CONFIG_PATH PINGVIN_INSTALL_SCHEMA_PATH)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/cmake/pingvin-config-version.cmake
    VERSION ${PINGVIN_VERSION_STRING}
    COMPATIBILITY AnyNewerVersion)

install(
    FILES
    ${CMAKE_BINARY_DIR}/cmake/pingvin-config.cmake
    ${CMAKE_BINARY_DIR}/cmake/pingvin-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pingvin)

# Create package
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
include(${CMAKE_SOURCE_DIR}/cmake/pingvin_cpack.cmake)
if(CPACK_GENERATOR)
    message(STATUS "Found CPack generators: ${CPACK_GENERATOR}")
    configure_file("${CMAKE_SOURCE_DIR}/cmake/cpack_options.cmake.in" ${PINGVIN_CPACK_CFG_FILE} @ONLY)
    set(CPACK_PROJECT_CONFIG_FILE ${PINGVIN_CPACK_CFG_FILE})
    include(CPack)
endif()
