ENABLE_TESTING()

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( ${GTEST_INCLUDE_DIRS}
                     ${CMAKE_BINARY_DIR}/apps/gadgetron
                     ${ACE_INCLUDE_DIR} 
                     ${Boost_INCLUDE_DIR}
                     ${FFTW3_INCLUDE_DIR}
                     ${ISMRMRD_INCLUDE_DIR}
                     ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
                     ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/math
                     ${CMAKE_SOURCE_DIR}/toolboxes/operators
                     ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/solvers
                     ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
                     ${CMAKE_SOURCE_DIR}/apps/gadgetron 
                     ${CMAKE_SOURCE_DIR}/toolboxes/core
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                     ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools 
                     ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow
                     ${CMAKE_SOURCE_DIR}/toolboxes/registration/optical_flow/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools
                     ${CMAKE_SOURCE_DIR}/toolboxes/gadgettools/ismrmrd )

link_libraries(optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY} 
                ${GTEST_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${ISMRMRD_LIBRARIES} 
                gadgetron_toolbox_cpucore 
                gadgetron_toolbox_cpucore_math 
                # cpureg
                gadgetron_toolbox_cpufft
                gadgetron_toolbox_gtplus 
                gadgetron_toolbox_gadgettools 
                gadgetron_toolbox_mri_core 
                gadgetronPlus)

if (CUDA_FOUND)
    link_libraries(gadgetron_toolbox_gtplus gadgetron_toolbox_gpuparallelmri gadgetron_toolbox_gpucore)
endif (CUDA_FOUND)

add_executable(gtplus_ut_grappa 
    gtplus_ut.cpp 
    grappa_test.cpp )

#add_test(gtplus_ut gtplus_ut_util)
