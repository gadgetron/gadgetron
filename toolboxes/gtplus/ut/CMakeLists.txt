ENABLE_TESTING()

if(WIN32)
    link_directories(${Boost_LIBRARY_DIRS})
endif(WIN32)

include_directories( ${GTEST_INCLUDE_DIRS} 
                     ${CMAKE_SOURCE_DIR}/gadgets/core 
                     ${ACE_INCLUDE_DIR} 
                     ${Boost_INCLUDE_DIR}
                     ${FFTW3_INCLUDE_DIR}
                     ${ISMRMRD_INCLUDE_DIR}
                     ${ISMRMRD_XSD_INCLUDE_DIR}
                     ${XSD_INCLUDE_DIR}
                     ${XERCESC_INCLUDE_DIR}
                     ${CMAKE_SOURCE_DIR}/dependencies/tinyxml
                     ${CMAKE_SOURCE_DIR}/toolboxes/mri/pmri/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/nfft/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/gpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/image
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/hostutils
                     ${CMAKE_SOURCE_DIR}/toolboxes/core/cpu/arma_math
                     ${CMAKE_SOURCE_DIR}/toolboxes/operators
                     ${CMAKE_SOURCE_DIR}/toolboxes/operators/cpu
                     ${CMAKE_SOURCE_DIR}/toolboxes/solvers
                     ${CMAKE_SOURCE_DIR}/toolboxes/solvers/cpu
                     ${CMAKE_SOURCE_DIR}/gadgets/core
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
                cpucore 
                cpucore_math 
                # cpureg
                cpufft
                gtplus 
                gadgettools 
                gadgetronPlus)

if (CUDA_FOUND AND CULA_FOUND)
    link_libraries(gtplus gpuparallelmri gpucore)
endif (CUDA_FOUND AND CULA_FOUND)

if (CUDA_FOUND AND CULA_FOUND)
    add_executable(gtplus_ut_grappa 
        gtplus_ut.cpp 
        grappa_test.cpp )
endif (CUDA_FOUND AND CULA_FOUND)

#add_test(gtplus_ut gtplus_ut_util)
